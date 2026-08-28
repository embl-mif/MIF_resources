"""Generic helpers for tile/ROI-based reading and writing of Zarr stores.
No tool-specific dependencies.
"""

import numpy as np


AXIS_NAMES = ("t", "c", "z", "y", "x")


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
