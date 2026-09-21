#!/usr/bin/env python3

"""Relabel one tile's segmentation with the global multicut mapping and
write its non-overlapping core region into the shared merged OME-Zarr store
(WRITE_RELABELED_TILE).

Runs once per tile, in parallel across every tile of the same source image,
each writing into a disjoint region of the store `PREALLOCATE_ZARR`
already created. A tile's own zarr array is already cropped back to its
base ROI by `CELLPOSE_SAM_ZARR` (halo is read-context only, see
`workflow_notes.md`), so this only needs to trim the *shared* band with
whichever neighbor(s) it has -- that band is always assigned to the
earlier tile along each axis (`--has-right`/`--has-bottom`/`--has-behind`
say whether this tile *is* an earlier tile on that axis), never split, so
no double-writes happen without needing any inter-tile coordination beyond
what's already in `tile_neighbors.csv`.

The write step reads the target store's own level-0 dataset path off its
on-disk attrs rather than assuming a fixed name (e.g. "s0") -- this store is
now always built by `preallocate_zarr.py --shape/--axes/--scale/--chunks
--translation` (mode 2, single-level; the old bespoke OMEZarrStore-based
preallocate_empty_zarr.py is retired), whose ngio-built naming mirrors
whatever axis/level naming the caller gave it and may nest multiscales
under an "ome" key for NGFF v0.5. Only the chunk-alignment helpers from
`ome_zarr_store.py` are reused here -- the `OMEZarrStore` class itself is
not, since its `.create()`-built layout is no longer the only supported
skeleton.
"""

import argparse
import json

import numpy as np
import zarr

from zarr_tile_utils import _read_tile_rows, compute_output_geometry
from label_matching import _load_full, _offset_labels, relabel_segmentation
from ome_zarr_store import check_grid_alignment, make_slice, storage_grid


def _multiscales_datasets(attrs) -> list[dict]:
    """Return the "datasets" list from a store's multiscales metadata,
    handling both the flat layout ome_zarr_store.py writes (multiscales at
    the attrs top level, both NGFF versions) and the NGFF v0.5-spec layout
    a real OME-Zarr library like ngio may write (multiscales nested under
    an "ome" key) -- same dual-handling precedent as
    `cellpose_sam_zarr.py`'s own `detect_ome_zarr_version`.
    """
    ome = attrs.get("ome")
    if isinstance(ome, dict) and "multiscales" in ome:
        multiscales = ome["multiscales"]
    else:
        multiscales = attrs.get("multiscales")

    if not multiscales:
        raise ValueError("Target store has no 'multiscales' metadata (checked both "
                          "top-level and 'ome'-nested locations).")

    return multiscales[0]["datasets"]


def _write_core_region(store_path: str, position: tuple[int, ...], data: np.ndarray) -> None:
    """Open the target store's level-0 array and write `data` at `position`,
    after checking the write lands on the array's own storage grid (see
    `check_grid_alignment` -- the boundary tile at each axis' true edge is
    exempted automatically, same rule `OMEZarrStore.write` enforced)."""
    root = zarr.open(store_path, mode="a")
    level0_path = _multiscales_datasets(root.attrs)[0]["path"]
    array = root[level0_path]

    shape = data.shape
    if not check_grid_alignment(position, shape, storage_grid(array), array.shape):
        raise ValueError(
            f"Tile write at position={position} shape={shape} is not aligned to the "
            f"storage grid of '{level0_path}' (grid={storage_grid(array)}, "
            f"array shape={array.shape})."
        )

    array[make_slice(position, shape)] = data


def _core_region(size: int, step: int, has_forward_neighbor: bool) -> tuple[slice, int]:
    """One axis' slice into this tile's own (local) array, and the resulting
    region length. The shared band with the next tile (width `overlap`) is
    dropped here -- the next tile will start its own core region right where
    this one ends, so nothing in between is ever written twice."""
    if has_forward_neighbor:
        return slice(0, step), step
    return slice(None), size


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument("--tile-label-path", type=str, required=True)
    parser.add_argument("--tiles-csv", type=str, required=True)
    parser.add_argument("--store-path", type=str, required=True)
    parser.add_argument("--global-mapping", type=str, required=True)

    parser.add_argument("--tile-index", type=int, required=True)
    parser.add_argument("--label-offset-step", type=int, required=True)

    parser.add_argument("--x-min", type=int, required=True)
    parser.add_argument("--y-min", type=int, required=True)
    parser.add_argument("--z-min", type=int, default=None)

    parser.add_argument("--tile-size-x", type=int, required=True)
    parser.add_argument("--tile-size-y", type=int, required=True)
    parser.add_argument("--tile-size-z", type=int, default=None)

    parser.add_argument("--tile-overlap-x", type=int, required=True)
    parser.add_argument("--tile-overlap-y", type=int, required=True)
    parser.add_argument("--tile-overlap-z", type=int, default=None)

    parser.add_argument("--has-right", action="store_true")
    parser.add_argument("--has-bottom", action="store_true")
    parser.add_argument("--has-behind", action="store_true")

    parser.add_argument("--background", type=int, default=0)
    parser.add_argument("--output", type=str, required=True)

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    seg = _load_full(args.tile_label_path)
    ndim = seg.ndim

    if ndim == 3:
        if args.z_min is None or args.tile_size_z is None or args.tile_overlap_z is None:
            raise ValueError("--z-min/--tile-size-z/--tile-overlap-z are required for a 3D tile")
        starts = (args.z_min, args.y_min, args.x_min)
        steps = (
            args.tile_size_z - args.tile_overlap_z,
            args.tile_size_y - args.tile_overlap_y,
            args.tile_size_x - args.tile_overlap_x,
        )
        has_forward = (args.has_behind, args.has_bottom, args.has_right)
    elif ndim == 2:
        starts = (args.y_min, args.x_min)
        steps = (
            args.tile_size_y - args.tile_overlap_y,
            args.tile_size_x - args.tile_overlap_x,
        )
        has_forward = (args.has_bottom, args.has_right)
    else:
        raise ValueError(f"Unsupported segmentation ndim={ndim}, expected 2 or 3")

    # tile_index * label_offset_step must match what Module 4 used for this
    # same tile when building the RAG, so this tile's own labels land on the
    # same node IDs `global_label_mapping.json`'s keys already use.
    offset_seg = _offset_labels(seg, args.tile_index * args.label_offset_step, args.background)

    with open(args.global_mapping, "r", encoding="utf-8") as handle:
        raw_mapping = json.load(handle)
    label_mapping = {int(key): int(value) for key, value in raw_mapping.items()}

    relabeled = relabel_segmentation(offset_seg, label_mapping, args.background)

    local_slices = []
    region_shape = []
    for size, step, forward in zip(seg.shape, steps, has_forward):
        local_slice, length = _core_region(size, step, forward)
        local_slices.append(local_slice)
        region_shape.append(length)

    core = relabeled[tuple(local_slices)]

    rows = _read_tile_rows(args.tiles_csv)
    _axis_order, _shape, _chunks, global_min = compute_output_geometry(rows)
    position = tuple(start - origin for start, origin in zip(starts, global_min))

    _write_core_region(args.store_path, position, core.astype(np.uint32))

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(args.store_path)

    print(f"Wrote tile at position={position} shape={core.shape} to {args.store_path}")


if __name__ == "__main__":
    main()
