#!/usr/bin/env python3

"""Compute the --start-level/--<axis>-min/--<axis>-max/--chunks flags that
let preallocate_zarr.py build the merged-segmentation OME-Zarr store as a
real crop of the original image's own pyramid.

Prints (to stdout, one line, space-separated) rather than building anything
itself, so it can be dropped into a Nextflow process's script block via
command substitution ahead of the actual preallocate_zarr.py call:

    params=$(calculate_zarr_params.py --tiles-csv "${tiles_csv}" \
        --source-zarr "${source_zarr}")
    preallocate_zarr.py --input-zarr "${source_zarr}" \
        --output-zarr "${output_zarr}" --drop-axes c \
        --dtype uint32 --data-type label $params
"""

import argparse

import zarr

from ome_zarr_store import OMEZarrMetadata
from zarr_tile_utils import _read_tile_rows, compute_output_geometry


def calculate_params(tiles_csv: str, source_zarr: str) -> dict:
    rows = _read_tile_rows(tiles_csv)
    if not rows:
        raise ValueError(f"{tiles_csv} has no tile rows")

    # resolution_level is a plain CSV column, constant across every row of one
    # source image (written once per image by tile_constructor.py) -- read it
    # here rather than threading it in separately as another CLI argument.
    resolution_level = int(rows[0]["resolution_level"])

    tile_axis_order, shape, chunks, global_min = compute_output_geometry(rows)
    tile_bounds = {
        axis: (global_min[i], global_min[i] + shape[i])
        for i, axis in enumerate(tile_axis_order)
    }
    tile_chunks = dict(zip(tile_axis_order, chunks))

    # The tile CSV only ever carries spatial (z/y/x) tile bounds -- tiles are
    # never split along other axes such as 't' (CELLPOSE_SAM_ZARR_MIF_TILED
    # loops over every 't' frame itself within a single tile task, and no two
    # tile tasks ever write overlapping 't' ranges of the same x/y/z region).
    # Any axis the source actually has that the tile CSV doesn't know about
    # must pass through uncropped, with one chunk spanning its whole extent --
    # there's no cross-task write boundary to protect on that axis.
    src_group = zarr.open(source_zarr, mode="r")
    src_axis_names = tuple(axis["name"] for axis in OMEZarrMetadata(src_group).axes)
    level_shape = src_group[str(resolution_level)].shape

    # Final axis order matches the SOURCE's own axis order (what
    # preallocate_zarr.py's output_axis_names will be), not the tile CSV's
    # z/y/x order, and excludes 'c' -- callers of this script always pass
    # --drop-axes c to preallocate_zarr.py downstream of these flags.
    axis_order = tuple(axis for axis in src_axis_names if axis != "c")
    bounds = {axis: tile_bounds[axis] for axis in axis_order if axis in tile_bounds}
    chunks_by_axis = {
        axis: (
            tile_chunks[axis]
            if axis in tile_chunks
            else level_shape[src_axis_names.index(axis)]
        )
        for axis in axis_order
    }

    return {
        "resolution_level": resolution_level,
        "axis_order": axis_order,
        "bounds": bounds,
        "chunks": chunks_by_axis,
    }


def format_as_flags(params: dict) -> str:
    """Render computed params as preallocate_zarr.py's own explicit-geometry
    flags -- comma-joined per flag, space-separated between flags, safe for
    a shell to word-split via `$(...)` command substitution."""
    flags = [f"--start-level {params['resolution_level']}"]
    for axis in params["axis_order"]:
        bounds = params["bounds"].get(axis)
        if bounds is not None:
            start, end = bounds
            flags.append(f"--{axis}-min {start}")
            flags.append(f"--{axis}-max {end}")
    flags.append(f"--chunks {','.join(str(params['chunks'][axis]) for axis in params['axis_order'])}")
    return " ".join(flags)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute the merged-segmentation OME-Zarr store's --input-zarr-mode "
                     "flags (--start-level / crop bounds / --chunks) from its tile CSV's "
                     "own ROI bounds and resolution_level, plus the source zarr's own "
                     "axes for anything the tile CSV doesn't cover (e.g. 't')."
    )

    parser.add_argument("--tiles-csv", type=str, required=True)
    parser.add_argument("--source-zarr", type=str, required=True,
                         help="The same store passed as preallocate_zarr.py's "
                              "--input-zarr -- read here to fill in axes (e.g. 't') "
                              "the tile CSV itself has no bounds for.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    params = calculate_params(args.tiles_csv, args.source_zarr)
    print(format_as_flags(params))


if __name__ == "__main__":
    main()
