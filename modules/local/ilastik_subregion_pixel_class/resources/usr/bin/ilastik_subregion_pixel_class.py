#!/usr/bin/env python3
import argparse
import fcntl
import json
import math
import os
import sys

import numpy as np
import zarr
import vigra
import h5py

from collections import OrderedDict
from pathlib import Path

os.environ.setdefault("LAZYFLOW_THREADS", "2")
os.environ.setdefault("LAZYFLOW_TOTAL_RAM_MB", "4000")

from ilastik import app
from ilastik.applets.dataSelection.opDataSelection import PreloadedArrayDatasetInfo


def _parse_optional_int(val):
    """Safely parse integer arguments that may be passed as empty strings from shell/Nextflow."""
    if val is None or val == "" or str(val).strip() == "":
        return None
    return int(val)

def extract_ilp_labels(ilp_path: str):
    """Extract label names directly from the Ilastik project file. Count determines class channels."""
    if not os.path.exists(ilp_path):
        raise FileNotFoundError(f"Ilastik project file not found: {ilp_path}")

    with h5py.File(ilp_path, "r") as f:
        if "PixelClassification/LabelNames" not in f:
            raise KeyError("Could not find '/PixelClassification/LabelNames' in the .ilp file.")

        raw_labels = f["PixelClassification/LabelNames"][()]
        return [l.decode("utf-8") if isinstance(l, bytes) else str(l) for l in raw_labels]


def _ilp_max_sigma(f: h5py.File) -> float:
    """Extract largest active sigma from selected feature scales."""
    if "FeatureSelections/Scales" not in f or "FeatureSelections/SelectionMatrix" not in f:
        return 1.6  # Default fallback sigma

    scales = np.asarray(f["FeatureSelections/Scales"][()], dtype=float).ravel()
    matrix = np.asarray(f["FeatureSelections/SelectionMatrix"][()])

    if matrix.size == 0 or not np.any(matrix):
        return 1.6

    active_indices = np.where(matrix.any(axis=0))[0]
    if len(active_indices) == 0:
        return 1.6

    return float(np.max(scales[active_indices]))


def inspect_ilp(project_path: str):
    """Inspect .ilp project to extract 2D/3D mode and calculate recommended halo."""
    with h5py.File(project_path, "r") as f:
        is_2d = False
        if "FeatureSelections/ComputeIn2d" in f:
            val = f["FeatureSelections/ComputeIn2d"][()]
            is_2d = bool(np.all(val) or np.any(val))

        max_sigma = _ilp_max_sigma(f)
        recommended_halo = int(math.ceil(3.5 * max_sigma))

    return is_2d, recommended_halo

def validate_grid_alignment(offsets, tile_shapes, zarr_array, axis_names=("t", "c", "z", "y", "x")):
    """
    Vectorized ND chunk alignment validation strictly for Zarr v2. Similar to https://github.com/embl-ic/squirrel/blob/main/squirrel/library/ome_zarr.py but raises descriptive ValueErrors identifying the exact failing axis.
    """
    offsets = np.asarray(offsets)
    tile_shapes = np.asarray(tile_shapes)
    chunks = np.asarray(zarr_array.chunks)
    full_shapes = np.asarray(zarr_array.shape)

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

def process_single_tile(tile_spec, in_store, out_store, shell, args, is_2d_model, halo):
    shape_t, shape_c, shape_z, shape_y, shape_x = in_store.shape

    # Extract tile coordinates
    raw_z_min = tile_spec.get("z_min", args.z_min)
    raw_z_max = tile_spec.get("z_max", args.z_max)
    
    is_3d_block = (raw_z_min is not None) and (raw_z_max is not None) and not is_2d_model

    y_min = tile_spec.get("y_min") if tile_spec.get("y_min") is not None else (args.y_min if args.y_min is not None else 0)
    y_max = tile_spec.get("y_max") if tile_spec.get("y_max") is not None else (args.y_max if args.y_max is not None else shape_y)
    x_min = tile_spec.get("x_min") if tile_spec.get("x_min") is not None else (args.x_min if args.x_min is not None else 0)
    x_max = tile_spec.get("x_max") if tile_spec.get("x_max") is not None else (args.x_max if args.x_max is not None else shape_x)

    z_min = raw_z_min if is_3d_block else (raw_z_min if raw_z_min is not None else 0)
    z_max = raw_z_max if is_3d_block else z_min + 1

    # Check boundaries against chunking configuration using vectorized validator
    validate_grid_alignment(
        offsets=(z_min, y_min, x_min),
        tile_shapes=(z_max - z_min, y_max - y_min, x_max - x_min),
        zarr_array=out_store,
        axis_names=("z", "y", "x")
    )

    if is_3d_block:
        z_halo = halo
        pz_min = max(0, z_min - z_halo)
        pz_max = min(shape_z, z_max + z_halo)
        z_read_slice = slice(pz_min, pz_max)

        cz_start = z_halo if z_min > 0 else 0
        cz_end = cz_start + (z_max - z_min)
        z_crop_slice = slice(cz_start, cz_end)
        z_write_slice = slice(z_min, z_max)
    else:
        z_read_slice = slice(z_min, z_max)
        z_crop_slice = slice(0, 1)
        z_write_slice = slice(z_min, z_max)

    py_min = max(0, y_min - halo)
    py_max = min(shape_y, y_max + halo)
    px_min = max(0, x_min - halo)
    px_max = min(shape_x, x_max + halo)

    t_index = tile_spec.get("t_index", args.t_index)

    # Construct strict 5D read slice for [t, c, z, y, x]
    in_read_slices = (
        slice(t_index, t_index + 1),  # t
        slice(0, 1),                  # c
        z_read_slice,                 # z
        slice(py_min, py_max),        # y
        slice(px_min, px_max)         # x
    )

    # Read tile data
    raw_chunk = in_store[tuple(in_read_slices)]

    # Run inference on preloaded shell
    tagged_chunk = vigra.taggedView(raw_chunk, "tczyx")
    role_data = OrderedDict([("Raw Data", [PreloadedArrayDatasetInfo(preloaded_array=tagged_chunk)])])

    predictions = shell.workflow.batchProcessingApplet.run_export(role_data, export_to_array=True)
    prob_chunk = predictions[0]

    # Resolve channel axis dynamically
    labels = extract_ilp_labels(args.project)
    n_classes = len(labels)

    expected_spatial = (
        z_read_slice.stop - z_read_slice.start,
        py_max - py_min,
        px_max - px_min,
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

    if prob_chunk.ndim == 5 and prob_chunk.shape[0] == 1:
        prob_chunk = prob_chunk[0]
        matching_c_axis -= 1

    if matching_c_axis != 0:
        prob_chunk = np.moveaxis(prob_chunk, matching_c_axis, 0)

    # Crop halos and write to target Zarr location
    cy_start = halo if y_min > 0 else 0
    cy_end = cy_start + (y_max - y_min)
    cx_start = halo if x_min > 0 else 0
    cx_end = cx_start + (x_max - x_min)

    crop_slices = (
        slice(None),
        z_crop_slice,
        slice(cy_start, cy_end),
        slice(cx_start, cx_end),
    )
    cropped_probs = prob_chunk[crop_slices]

    if out_store.dtype == np.uint8 and cropped_probs.dtype != np.uint8:
        cropped_probs = np.rint(np.clip(cropped_probs, 0, 1) * 255).astype(np.uint8)

    write_slices = (
        slice(t_index, t_index + 1),
        slice(None),
        z_write_slice,
        slice(y_min, y_max),
        slice(x_min, x_max),
    )
    out_store[write_slices] = cropped_probs[np.newaxis, ...]
    print(f"[Worker] Finished Tile Y:[{y_min}:{y_max}] X:[{x_min}:{x_max}] Z:[{z_min}:{z_max}]")


def ensure_preallocated(input_path: str, ilp_path: str, output_path: str, lock_path: str):
    """
    Double-checked locking pattern using Linux standard library fcntl.
    Ensures the skeleton Zarr store is created exactly once, even if hundreds
    of worker nodes start at the same second.
    """
    marker_file = os.path.join(output_path, ".zgroup")

    # 1st Check: Fast path if store already exists (bypasses locking overhead)
    if os.path.exists(marker_file):
        return

    # Create lock file directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)

    print(f"[Prealloc Lock] Waiting for exclusive lock on {lock_path}...")
    with open(lock_path, "w") as lock_file_obj:
        # Block until acquiring exclusive lock across cluster filesystem
        fcntl.flock(lock_file_obj, fcntl.LOCK_EX)
        try:
            # 2nd Check: Verify another node didn't build it while this job was waiting
            if os.path.exists(marker_file):
                print("[Prealloc Lock] Store created by another worker while waiting. Proceeding.")
                return

            print("[Prealloc Lock] Acquired lock. Preallocating empty OME-Zarr skeleton...")

            labels = extract_ilp_labels(ilp_path)
            num_channels = len(labels)
            print(f"[Prealloc Lock] Detected {num_channels} target classes from .ilp: {labels}")

            src_group = zarr.open(input_path, mode="r")
            store = zarr.DirectoryStore(output_path)
            dest_group = zarr.group(store=store, overwrite=True)

            zattrs = dict(src_group.attrs)
            multiscales = zattrs.get("multiscales", [])
            if not multiscales:
                raise ValueError(f"Input Zarr '{input_path}' lacks 'multiscales' metadata.")

            axes = multiscales[0].get("axes", [])
            channel_axis_index = next(
                (idx for idx, axis in enumerate(axes) if (axis.get("name") if isinstance(axis, dict) else axis) == "c"),
                1
            )

            # Update OMERO visualizer channel colors
            if "omero" in zattrs:
                palette = ["FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "FFFFFF"]
                zattrs["omero"]["channels"] = [
                    {
                        "active": True,
                        "coefficient": 1.0,
                        "color": palette[i % len(palette)],
                        "family": "linear",
                        "inverted": False,
                        "label": label_name,
                        "window": {"end": 1.0, "max": 1.0, "min": 0.0, "start": 0.0},
                    }
                    for i, label_name in enumerate(labels)
                ]

            dest_group.attrs.update(zattrs)

            # Create empty arrays for resolution levels
            for dataset_meta in multiscales[0]["datasets"]:
                path_name = dataset_meta["path"]
                src_array = src_group[path_name]

                new_shape = list(src_array.shape)
                new_shape[channel_axis_index] = num_channels

                # Match chunk size for channels to num_channels for atomic write operations
                new_chunks = list(src_array.chunks)
                new_chunks[channel_axis_index] = num_channels

                print(f"[Prealloc Lock] Creating scale '{path_name}': shape {src_array.shape} -> {new_shape}")

                dest_group.create_dataset(
                    name=path_name,
                    shape=new_shape,
                    chunks=new_chunks,
                    dtype="float32",
                    compressor=src_array.compressor,
                    fill_value=0.0,
                )

            print("[Prealloc Lock] Preallocation complete. Releasing lock.")

        finally:
            # Guarantee lock release
            fcntl.flock(lock_file_obj, fcntl.LOCK_UN)
def main():
    parser = argparse.ArgumentParser(description="Ilastik tile worker with built-in native file lock preallocation.")
    parser.add_argument("--input-zarr", required=True, help="Path to input raw Zarr store")
    parser.add_argument("--output-zarr", required=True, help="Path to output probability Zarr store")
    parser.add_argument("--project", required=True, help="Path to Ilastik project file (.ilp)")

    parser.add_argument("--t-index", type=int, default=0, help="Time index")
    parser.add_argument("--z-min", type=_parse_optional_int, default=None, help="Min Z index")
    parser.add_argument("--z-max", type=_parse_optional_int, default=None, help="Max Z index")

    # Made spatial boundaries optional to allow full-image processing when unsupplied
    parser.add_argument("--y-min", type=_parse_optional_int, default=None, help="Min Y index (defaults to 0)")
    parser.add_argument("--y-max", type=_parse_optional_int, default=None, help="Max Y index (defaults to full height)")
    parser.add_argument("--x-min", type=_parse_optional_int, default=None, help="Min X index (defaults to 0)")
    parser.add_argument("--x-max", type=_parse_optional_int, default=None, help="Max X index (defaults to full width)")

    parser.add_argument("--halo", type=_parse_optional_int, default=None, help="Halo pixel size")
    parser.add_argument("--done-file", type=str, default=None, help="Output JSON path upon completion")

    args = parser.parse_args()
    
    # lock-protected preallocation
    lock_file = f"{args.output_zarr}.lock"
    ensure_preallocated(args.input_zarr, args.project, args.output_zarr, lock_file)

    # Open Zarr stores
    in_root = zarr.open(args.input_zarr, mode="r")
    in_store = in_root["0"] if isinstance(in_root, zarr.hierarchy.Group) and "0" in in_root else in_root

    out_root = zarr.open(args.output_zarr, mode="r+")
    out_store = out_root["0"] if isinstance(out_root, zarr.hierarchy.Group) and "0" in out_root else out_root

    # Determine if project is 2D or 3D
    is_2d_model, auto_halo = inspect_ilp(args.project)
    halo = args.halo if args.halo is not None else auto_halo

    print("[Worker] Initializing Ilastik Headless Shell...")
    ilastik_args = app.parse_args([])
    ilastik_args.headless = True
    ilastik_args.project = args.project
    shell = app.main(ilastik_args)

    # Figure out tiles to process
    if args.tiles_json:
        with open(args.tiles_json, "r") as f:
            tiles = json.load(f)
    else:
        # Fallback to single tile
        tiles = [{
            "y_min": args.y_min, "y_max": args.y_max,
            "x_min": args.x_min, "x_max": args.x_max,
            "z_min": args.z_min, "z_max": args.z_max,
            "t_index": args.t_index
        }]

    for tile in tiles:
        process_single_tile(tile, in_store, out_store, shell, args, is_2d_model, halo)

    # Write single done file for the batch
    if args.done_file:
        done_path = Path(args.done_file)
        done_path.parent.mkdir(parents=True, exist_ok=True)
        with open(done_path, "w") as f:
            json.dump({"processed_count": len(tiles), "status": "COMPLETED"}, f, indent=2)

if __name__ == "__main__":
    main()
