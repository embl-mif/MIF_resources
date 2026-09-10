"""Generic helpers for tile/ROI-based reading and writing of Zarr stores.

compute_output_geometry/_read_tile_rows/_read_input_calibration moved here
from the now-retired preallocate_empty_zarr.py: the stitched multi-tile
merge store is built via the general-purpose preallocate_zarr.py
--shape/--axes/--scale/--chunks/--translation mode instead of a
bespoke OMEZarrStore-based script, and this is the geometry math that feeds
that mode. write_relabeled_tile.py also needs compute_output_geometry/
_read_tile_rows for its own write-position math, so this is the one place
both the merge-store-creation step and the per-tile write step can import
it from, without either depending on preallocate_zarr.py itself.
"""

import csv

import numpy as np
import zarr

from ome_zarr_store import detect_ome_version, detect_zarr_format


AXIS_NAMES = ("t", "c", "z", "y", "x")


def _read_tile_rows(tiles_csv: str) -> list[dict[str, str]]:
    with open(tiles_csv, "r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _has_z(rows: list[dict[str, str]]) -> bool:
    return bool(rows[0].get("z_min", "").strip())


def _bounds(rows: list[dict[str, str]], axis: str) -> tuple[int, int]:
    mins = [int(row[f"{axis}_min"]) for row in rows]
    maxs = [int(row[f"{axis}_max"]) for row in rows]
    return min(mins), max(maxs)


def compute_output_geometry(
    rows: list[dict[str, str]],
) -> tuple[list[str], tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    """Derive the merged output's geometry from every tile row of one source image.

    Returns
    -------
    axis_order : list[str]
        ``["z", "y", "x"]`` or ``["y", "x"]``, matching every other tuple below.
    shape : tuple[int, ...]
        Union of every tile's ROI bounds -- the output array's own shape.
    chunks : tuple[int, ...]
        Tiling stride (`tile_size - tile_overlap`) per axis -- the output
        store's chunk size.
    global_min : tuple[int, ...]
        The output array's origin, in the original image's own voxel
        coordinates. Every tile's write position is its own ROI start minus
        this.
    """
    axis_order = ["z", "y", "x"] if _has_z(rows) else ["y", "x"]

    bounds = {axis: _bounds(rows, axis) for axis in axis_order}
    shape = tuple(bounds[axis][1] - bounds[axis][0] for axis in axis_order)
    global_min = tuple(bounds[axis][0] for axis in axis_order)

    first = rows[0]
    chunks = tuple(
        int(first[f"tile_size_{axis}"]) - int(first[f"tile_overlap_{axis}"]) for axis in axis_order
    )

    return axis_order, shape, chunks, global_min


def _read_input_calibration(
    input_image_path: str, resolution_level: int, axis_order: list[str]
) -> tuple[str, int, tuple[float, ...]]:
    """Detect the input image's OME version/Zarr format, and read its
    per-axis voxel scale at `resolution_level` -- for calibrating a merge
    output store whose shape comes from `compute_output_geometry`, not from
    this image."""
    root = zarr.open(input_image_path, mode="r")
    ome_version = detect_ome_version(root)
    zarr_format = detect_zarr_format(root)

    ms = root.attrs["multiscales"][0]
    dataset = ms["datasets"][resolution_level]
    axes = [axis["name"] for axis in ms["axes"]]

    scale_by_axis = {"x": 1.0, "y": 1.0, "z": 1.0}
    for transform in dataset["coordinateTransformations"]:
        if transform.get("type") == "scale":
            for name, value in zip(axes, transform["scale"]):
                if name in scale_by_axis:
                    scale_by_axis[name] = float(value)

    resolution = tuple(scale_by_axis[axis] for axis in axis_order)
    return ome_version, zarr_format, resolution


def parse_optional_int(val):
    """turns empty string or None into None, otherwise int. handles nextflow passing empty strings for unset args."""
    if val is None or val == "" or str(val).strip() == "":
        return None
    return int(val)


def expand_with_halo(start, end, halo_voxels, axis_size):
    """
    Expand a voxel-space [start, end) region by a fixed halo margin on each
    side, clipped to the array's real extent. Returns (load_start, load_end).
    """
    if halo_voxels <= 0:
        return start, end
    load_start = max(0, start - halo_voxels)
    load_end = min(axis_size, end + halo_voxels)
    return load_start, load_end


def crop_slice_for_halo(base_start, base_end, load_start):
    """
    Given a halo-expanded read window starting at load_start, compute the
    slice (relative to that window) that recovers exactly [base_start, base_end).
    """
    crop_start = base_start - load_start
    crop_end = crop_start + (base_end - base_start)
    return slice(crop_start, crop_end)


def expand_and_crop_halo(bounds, halo_per_axis, shapes):
    """
    Given base region bounds, a per-axis halo margin, and the array's real
    shape, compute the halo-expanded read slice and the crop slice that
    recovers the original bounds from it -- for any number of axes at once.

    bounds:        dict {axis_name: (start, end)}
    halo_per_axis: dict {axis_name: halo_voxels}
    shapes:        dict {axis_name: axis_size}

    Returns (read_slices, crop_slices): dicts of slice objects, keyed the
    same as bounds. read_slices index the full array; crop_slices index
    the block that read_slices produced.
    """
    read_slices = {}
    crop_slices = {}
    for axis, (start, end) in bounds.items():
        halo_voxels = halo_per_axis.get(axis, 0)
        load_start, load_end = expand_with_halo(start, end, halo_voxels, shapes[axis])
        read_slices[axis] = slice(load_start, load_end)
        crop_slices[axis] = crop_slice_for_halo(start, end, load_start)
    return read_slices, crop_slices


def validate_grid_alignment(offsets, tile_shapes, zarr_array, axis_names=("t", "c", "z", "y", "x")):
    """
    Check a tile's offset/size line up with the store's chunk grid. Only compatible with Zarr v2.
    Same logic as check_grid_alignment https://github.com/embl-ic/squirrel/blob/main/squirrel/library/ome_zarr.py
    """
    offsets = np.asarray(offsets)
    tile_shapes = np.asarray(tile_shapes)

    axis_index = [AXIS_NAMES.index(n) for n in axis_names]
    chunks = np.asarray(zarr_array.chunks)[axis_index]
    full_shapes = np.asarray(zarr_array.shape)[axis_index]

    for p, s, g, ds, name in zip(offsets, tile_shapes, chunks, full_shapes, axis_names):
        if p % g != 0:
            raise ValueError(
                f"Chunk Alignment Conflict ('{name}'): min offset ({p}) "
                f"is not a multiple of Zarr v2 chunk size ({g})."
            )
        if (s % g != 0) and (p + s < ds):
            raise ValueError(
                f"Chunk Alignment Conflict ('{name}'): tile size ({s}) "
                f"is not a multiple of Zarr v2 chunk size ({g}) and does not end at dataset boundary."
            )
