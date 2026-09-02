#!/usr/bin/env python3
"""Ilastik headless tile worker. Reads an already-preallocated OME-Zarr
probability store (see PREALLOCATE_ZARR / EXTRACT_ILASTIK_CLASSES modules) and
writes ilastik predictions into it, tile by tile, over a given timepoint range.
"""

import argparse

import numpy as np
import zarr
import vigra

from ilastik import app
from ilastik.applets.dataSelection.opDataSelection import PreloadedArrayDatasetInfo

from ilp_utils import extract_ilp_classes, inspect_ilp
from zarr_tile_utils import parse_optional_int, validate_grid_alignment, expand_and_crop_halo


def process_single_tile(t_index,in_store, out_store, shell, args, is_2d_model, halo, n_classes):
    """read one tile from in_store, run ilastik on it, crop off the halo, write into out_store."""
    if in_store.ndim != 5:
        raise ValueError(
            f"This worker expects a 5D (t, c, z, y, x) input array, got shape {in_store.shape}."
        )
    shape_t, shape_c, shape_z, shape_y, shape_x = in_store.shape

    y_min = args.y_min if args.y_min is not None else 0
    y_max = args.y_max if args.y_max is not None else shape_y
    x_min = args.x_min if args.x_min is not None else 0
    x_max = args.x_max if args.x_max is not None else shape_x
    z_min = args.z_min if args.z_min is not None else 0
    z_max = args.z_max if args.z_max is not None else shape_z
    t_index = int(t_index)

    if not (0 <= z_min < z_max <= shape_z):
        raise ValueError(f"Bad z range [{z_min}, {z_max}) for extent {shape_z}.")
    if not (0 <= y_min < y_max <= shape_y):
        raise ValueError(f"Bad y range [{y_min}, {y_max}) for extent {shape_y}.")
    if not (0 <= x_min < x_max <= shape_x):
        raise ValueError(f"Bad x range [{x_min}, {x_max}) for extent {shape_x}.")

    validate_grid_alignment(
        offsets=(z_min, y_min, x_min),
        tile_shapes=(z_max - z_min, y_max - y_min, x_max - x_min),
        zarr_array=out_store,
        axis_names=("z", "y", "x")
    )

    z_halo_voxels = 0 if is_2d_model else halo

    read, crop = expand_and_crop_halo(
        bounds={"z": (z_min, z_max), "y": (y_min, y_max), "x": (x_min, x_max)},
        halo_per_axis={"z": z_halo_voxels, "y": halo, "x": halo},
        shapes={"z": shape_z, "y": shape_y, "x": shape_x},
    )

    in_read_slices = (
        slice(t_index, t_index + 1),
        slice(args.c_index, args.c_index + 1),
        read["z"], read["y"], read["x"],
    )

    raw_chunk = np.ascontiguousarray(in_store[in_read_slices])

    tagged_chunk = vigra.taggedView(raw_chunk, "tczyx")
    role_data = [{"Raw Data": PreloadedArrayDatasetInfo(preloaded_array=tagged_chunk)}]
    predictions = shell.workflow.batchProcessingApplet.run_export(role_data, export_to_array=True)
    prob_chunk = np.asarray(predictions[0])

    expected_spatial = (
        read["z"].stop - read["z"].start,
        read["y"].stop - read["y"].start,
        read["x"].stop - read["x"].start,
    )

    prob_shape = list(prob_chunk.shape)
    cand_channel_axes = [i for i, s in enumerate(prob_shape) if s == n_classes]

    matching_c_axis = None
    for c_idx in cand_channel_axes:
        rem_shape = tuple(s for i, s in enumerate(prob_shape) if i != c_idx)
        if len(rem_shape) == 4 and rem_shape[0] == 1:
            rem_shape = rem_shape[1:]

        if rem_shape == expected_spatial:
            matching_c_axis = c_idx
            break

    if matching_c_axis is None:
        raise ValueError(
            f"Unable to safely identify channel axis in Ilastik export shape {prob_chunk.shape}. "
            f"Expected spatial shape {expected_spatial} for {n_classes} classes."
        )

    if prob_chunk.ndim == 5 and prob_chunk.shape[0] == 1 and matching_c_axis != 0:
        prob_chunk = prob_chunk[0]
        matching_c_axis -= 1

    if matching_c_axis != 0:
        prob_chunk = np.moveaxis(prob_chunk, matching_c_axis, 0)

    if prob_chunk.ndim != 4:
        raise ValueError(
            f"Expected a 4D (c, z, y, x) prediction after reordering, got {prob_chunk.shape}."
        )

    cropped_probs = prob_chunk[(slice(None), crop["z"], crop["y"], crop["x"])]

    expected_write = (n_classes, z_max - z_min, y_max - y_min, x_max - x_min)
    if cropped_probs.shape != expected_write:
        raise ValueError(
            f"Cropped prediction has shape {cropped_probs.shape}, expected {expected_write}."
        )

    if out_store.dtype == np.uint8 and cropped_probs.dtype != np.uint8:
        cropped_probs = np.rint(np.clip(cropped_probs, 0, 1) * 255).astype(np.uint8)

    write_slices = (
        slice(t_index, t_index + 1),
        slice(None),
        slice(z_min, z_max),
        slice(y_min, y_max),
        slice(x_min, x_max),
    )
    out_store[write_slices] = cropped_probs[np.newaxis, ...]
    print(f"[Worker] Finished Tile T:{t_index} Z:[{z_min}:{z_max}] Y:[{y_min}:{y_max}] X:[{x_min}:{x_max}]")


def main():
    """cli entrypoint: run ilastik tile by tile over the given time range, against an
    already-preallocated output store (see PREALLOCATE_ZARR)."""
    parser = argparse.ArgumentParser(description="Ilastik headless tile worker.")
    parser.add_argument("--input-zarr", required=True, help="Path to input raw Zarr store")
    parser.add_argument("--output-zarr", required=True,
                         help="Path to the output probability Zarr store. Must already exist "
                              "-- see PREALLOCATE_ZARR.")
    parser.add_argument("--project", required=True, help="Path to Ilastik project file (.ilp)")
    parser.add_argument("--t-min", type=parse_optional_int, default=None,
                         help="First timepoint index to process (inclusive). Defaults to 0.")
    parser.add_argument("--t-max", type=parse_optional_int, default=None,
                         help="Last timepoint index to process (exclusive). Defaults to full extent.")
    parser.add_argument("--c-index", type=int, default=0, help="Input channel index fed to Ilastik")
    parser.add_argument("--z-min", type=parse_optional_int, default=None, help="Min Z index")
    parser.add_argument("--z-max", type=parse_optional_int, default=None, help="Max Z index")
    parser.add_argument("--y-min", type=parse_optional_int, default=None, help="Min Y index (defaults to 0)")
    parser.add_argument("--y-max", type=parse_optional_int, default=None, help="Max Y index (defaults to full height)")
    parser.add_argument("--x-min", type=parse_optional_int, default=None, help="Min X index (defaults to 0)")
    parser.add_argument("--x-max", type=parse_optional_int, default=None, help="Max X index (defaults to full width)")
    parser.add_argument("--halo", type=parse_optional_int, default=None, help="Halo pixel size")

    args = parser.parse_args()

    in_root = zarr.open(args.input_zarr, mode="r")
    in_store = in_root["0"] if isinstance(in_root, zarr.hierarchy.Group) and "0" in in_root else in_root

    out_root = zarr.open(args.output_zarr, mode="r+")
    out_store = out_root["0"] if isinstance(out_root, zarr.hierarchy.Group) and "0" in out_root else out_root

    is_2d_model, auto_halo = inspect_ilp(args.project)
    halo = args.halo if args.halo is not None else auto_halo

    n_classes = len(extract_ilp_classes(args.project))
    if out_store.shape[1] != n_classes:
        raise ValueError(
            f"Output store has {out_store.shape[1]} channels but the project has "
            f"{n_classes} classes -- it was built from a different .ilp."
        )
    print(f"[Worker] {n_classes} classes, halo={halo}, 2D features={is_2d_model}")

    print("[Worker] Initializing Ilastik Headless Shell...")
    ilastik_args = app.parse_args([])
    ilastik_args.headless = True
    ilastik_args.project = args.project
    shell = app.main(ilastik_args)

    shape_t = in_store.shape[0]
    t_start = args.t_min if args.t_min is not None else 0
    t_end = args.t_max if args.t_max is not None else shape_t

    if not (0 <= t_start < t_end <= shape_t):
        raise ValueError(f"Bad t range [{t_start}, {t_end}) for extent {shape_t}.")

    for t_index in range(t_start, t_end):
        process_single_tile(t_index, in_store, out_store, shell, args, is_2d_model, halo, n_classes)
    print(f"[Worker] Finished all {t_end - t_start} tile(s).")

if __name__ == "__main__":
    main()
