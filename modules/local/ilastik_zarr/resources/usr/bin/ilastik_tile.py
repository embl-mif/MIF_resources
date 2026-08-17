import argparse
import fcntl
import json
import math
import os
import shutil

import numpy as np
import zarr
import vigra
import h5py

from pathlib import Path

os.environ.setdefault("LAZYFLOW_THREADS", "2")
os.environ.setdefault("LAZYFLOW_TOTAL_RAM_MB", "4000")

from ilastik import app
from ilastik.applets.dataSelection.opDataSelection import PreloadedArrayDatasetInfo


AXIS_NAMES = ("t", "c", "z", "y", "x")


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
    """Inspect .ilp project to extract 2D/3D mode and calculate recommended halo from sigmas."""
    with h5py.File(project_path, "r") as f:
        is_2d = False
        if "FeatureSelections/ComputeIn2d" in f:
            val = f["FeatureSelections/ComputeIn2d"][()]
            is_2d = bool(np.any(val))

        max_sigma = _ilp_max_sigma(f)
        recommended_halo = int(math.ceil(3.5 * max_sigma))

    return is_2d, recommended_halo


def validate_grid_alignment(offsets, tile_shapes, zarr_array, axis_names=("t", "c", "z", "y", "x")):
    """
    Vectorized ND chunk alignment validation strictly for Zarr v2.
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


def process_single_tile(tile_spec, in_store, out_store, shell, args, is_2d_model, halo, n_classes):
     """Read one tile from in_store, run ilastik on it, crop off the halo, write into out_store."""
    if in_store.ndim != 5:
        raise ValueError(
            f"This worker expects a 5D (t, c, z, y, x) input array, got shape {in_store.shape}."
        )
    shape_t, shape_c, shape_z, shape_y, shape_x = in_store.shape

    def _pick(key, cli_value, default):
        value = tile_spec.get(key)
        if value is None:
            value = cli_value
        return default if value is None else int(value)

    y_min = _pick("y_min", args.y_min, 0)
    y_max = _pick("y_max", args.y_max, shape_y)
    x_min = _pick("x_min", args.x_min, 0)
    x_max = _pick("x_max", args.x_max, shape_x)
    z_min = _pick("z_min", args.z_min, 0)
    z_max = _pick("z_max", args.z_max, shape_z)
    t_index = int(tile_spec["t_index"])

    if not (0 <= z_min < z_max <= shape_z):
        raise ValueError(f"Bad z range [{z_min}, {z_max}) for extent {shape_z}.")
    if not (0 <= y_min < y_max <= shape_y):
        raise ValueError(f"Bad y range [{y_min}, {y_max}) for extent {shape_y}.")
    if not (0 <= x_min < x_max <= shape_x):
        raise ValueError(f"Bad x range [{x_min}, {x_max}) for extent {shape_x}.")

    # Check boundaries against chunking configuration using vectorized validator
    validate_grid_alignment(
        offsets=(z_min, y_min, x_min),
        tile_shapes=(z_max - z_min, y_max - y_min, x_max - x_min),
        zarr_array=out_store,
        axis_names=("z", "y", "x")
    )

    z_halo = 0 if is_2d_model else halo
    pz_min = max(0, z_min - z_halo)
    pz_max = min(shape_z, z_max + z_halo)
    py_min = max(0, y_min - halo)
    py_max = min(shape_y, y_max + halo)
    px_min = max(0, x_min - halo)
    px_max = min(shape_x, x_max + halo)

    # Crop offsets are derived from the actual clipped read window rather
    # than assuming a full halo was available on the low side.
    z_crop_slice = slice(z_min - pz_min, (z_min - pz_min) + (z_max - z_min))
    y_crop_slice = slice(y_min - py_min, (y_min - py_min) + (y_max - y_min))
    x_crop_slice = slice(x_min - px_min, (x_min - px_min) + (x_max - x_min))

    # Construct strict 5D read slice for [t, c, z, y, x]
    in_read_slices = (
        slice(t_index, t_index + 1),          # t
        slice(args.c_index, args.c_index + 1),  # c
        slice(pz_min, pz_max),                # z
        slice(py_min, py_max),                # y
        slice(px_min, px_max),                # x
    )

    # Read tile data
    raw_chunk = np.ascontiguousarray(in_store[in_read_slices])

    # Run inference on preloaded shell
    tagged_chunk = vigra.taggedView(raw_chunk, "tczyx")
    role_data = [{"Raw Data": PreloadedArrayDatasetInfo(preloaded_array=tagged_chunk)}]
    predictions = shell.workflow.batchProcessingApplet.run_export(role_data, export_to_array=True)
    prob_chunk = np.asarray(predictions[0])

    # Resolve channel axis dynamically
    expected_spatial = (
        pz_max - pz_min,
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

    if prob_chunk.ndim == 5 and prob_chunk.shape[0] == 1 and matching_c_axis != 0:
        prob_chunk = prob_chunk[0]
        matching_c_axis -= 1

    if matching_c_axis != 0:
        prob_chunk = np.moveaxis(prob_chunk, matching_c_axis, 0)

    if prob_chunk.ndim != 4:
        raise ValueError(
            f"Expected a 4D (c, z, y, x) prediction after reordering, got {prob_chunk.shape}."
        )

    # Crop halos and write to target Zarr location
    cropped_probs = prob_chunk[(slice(None), z_crop_slice, y_crop_slice, x_crop_slice)]

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


def _build_skeleton(input_path: str, ilp_path: str, output_path: str, dtype: str):
    """Create the empty OME-Zarr skeleton at output_path (which must not exist)."""
    labels = extract_ilp_labels(ilp_path)
    num_channels = len(labels)
    print(f"[Prealloc] Detected {num_channels} target classes from .ilp: {labels}")

    src_group = zarr.open(input_path, mode="r")
    dest_group = zarr.group(store=zarr.DirectoryStore(output_path), overwrite=True)

    zattrs = dict(src_group.attrs)
    multiscales = zattrs.get("multiscales", [])
    if not multiscales:
        raise ValueError(f"Input Zarr '{input_path}' lacks 'multiscales' metadata.")

    axes = multiscales[0].get("axes", [])
    channel_axis_index = next(
        (idx for idx, axis in enumerate(axes)
         if (axis.get("name") if isinstance(axis, dict) else axis) == "c"),
        None
    )

    if channel_axis_index is None:
        raise ValueError(
            f"Input axes {axes} contain no 'c' axis; cannot place {num_channels} "
            "probability channels."
        )

    max_val = 255.0 if dtype == "uint8" else 1.0
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
                "window": {"end": max_val, "max": max_val, "min": 0.0, "start": 0.0},
            }
            for i, label_name in enumerate(labels)
        ]

    dest_group.attrs.update(zattrs)

    for dataset_meta in multiscales[0]["datasets"]:
        path_name = dataset_meta["path"]
        src_array = src_group[path_name]

        new_shape = list(src_array.shape)
        new_shape[channel_axis_index] = num_channels

        # Match chunk size for channels to num_channels for atomic write operations
        new_chunks = list(src_array.chunks)
        new_chunks[channel_axis_index] = num_channels

        separator = getattr(src_array, "_dimension_separator", None) or "."

        print(f"[Prealloc] Creating scale '{path_name}': shape {src_array.shape} -> {new_shape}")

        dest_group.create_dataset(
            name=path_name,
            shape=new_shape,
            chunks=new_chunks,
            dtype=dtype,
            compressor=src_array.compressor,
            fill_value=0,
            dimension_separator=separator,
        )


def ensure_preallocated(input_path: str, ilp_path: str, output_path: str,
                        lock_path: str, dtype: str = "float32"):
    """
    Build the output skeleton once, safely, even with many jobs starting at the same time.

    Checks a marker file that only gets written once the skeleton is fully built, not just 
    when the store's created). So a job that shows up mid-build waits
    on the lock instead of grabbing a half-finished store, and a crashed build gets
    redone next time instead of getting stuck forever.
    """
    marker_file = f"{output_path}.complete"

    # 1st Check: fast path if the store is known-complete (no locking overhead)
    if os.path.exists(marker_file):
        return

    os.makedirs(os.path.dirname(os.path.abspath(lock_path)) or ".", exist_ok=True)

    print(f"[Prealloc Lock] Waiting for exclusive lock on {lock_path}...")
    with open(lock_path, "w") as lock_file_obj:
        # Block until acquiring exclusive lock across cluster filesystem
        fcntl.flock(lock_file_obj, fcntl.LOCK_EX)
        try:
            # 2nd Check: verify another node didn't build it while this job waited
            if os.path.exists(marker_file):
                print("[Prealloc Lock] Store created by another worker while waiting. Proceeding.")
                return

            print("[Prealloc Lock] Acquired lock. Preallocating empty OME-Zarr skeleton...")

            # Any store present without the marker is from an interrupted build.
            if os.path.exists(output_path):
                print("[Prealloc Lock] Removing incomplete store from a previous attempt.")
                shutil.rmtree(output_path)

            _build_skeleton(input_path, ilp_path, output_path, dtype)

            with open(marker_file, "w") as fh:
                json.dump({"status": "COMPLETE", "dtype": dtype}, fh)

            print("[Prealloc Lock] Preallocation complete. Releasing lock.")

        finally:
            # Guarantee lock release
            fcntl.flock(lock_file_obj, fcntl.LOCK_UN)


def main():
    """cli entrypoint: preallocate the output store if needed, then run ilastik tile by tile over the given time range."""
    parser = argparse.ArgumentParser(description="Ilastik tile worker with built-in native file lock preallocation.")
    parser.add_argument("--input-zarr", required=True, help="Path to input raw Zarr store")
    parser.add_argument("--output-zarr", required=True, help="Path to output probability Zarr store")
    parser.add_argument("--project", required=True, help="Path to Ilastik project file (.ilp)")
    parser.add_argument("--t-min", type=_parse_optional_int, default=None,
                         help="First timepoint index to process (inclusive). Defaults to 0.")
    parser.add_argument("--t-max", type=_parse_optional_int, default=None,
                         help="Last timepoint index to process (exclusive). Defaults to full extent.")
    parser.add_argument("--c-index", type=int, default=0, help="Input channel index fed to Ilastik")
    parser.add_argument("--z-min", type=_parse_optional_int, default=None, help="Min Z index")
    parser.add_argument("--z-max", type=_parse_optional_int, default=None, help="Max Z index")

    # Made spatial boundaries optional to allow full-image processing when unsupplied
    parser.add_argument("--y-min", type=_parse_optional_int, default=None, help="Min Y index (defaults to 0)")
    parser.add_argument("--y-max", type=_parse_optional_int, default=None, help="Max Y index (defaults to full height)")
    parser.add_argument("--x-min", type=_parse_optional_int, default=None, help="Min X index (defaults to 0)")
    parser.add_argument("--x-max", type=_parse_optional_int, default=None, help="Max X index (defaults to full width)")

    parser.add_argument("--halo", type=_parse_optional_int, default=None, help="Halo pixel size")
    parser.add_argument("--dtype", choices=["float32", "uint8"], default="float32",
                        help="Probability dtype for the preallocated store")
    parser.add_argument("--done-file", type=str, default=None, help="Output JSON path upon completion")

    args = parser.parse_args()

    # lock-protected preallocation
    lock_file = f"{args.output_zarr}.lock"
    ensure_preallocated(args.input_zarr, args.project, args.output_zarr, lock_file, args.dtype)

    # Open Zarr stores
    in_root = zarr.open(args.input_zarr, mode="r")
    in_store = in_root["0"] if isinstance(in_root, zarr.hierarchy.Group) and "0" in in_root else in_root

    out_root = zarr.open(args.output_zarr, mode="r+")
    out_store = out_root["0"] if isinstance(out_root, zarr.hierarchy.Group) and "0" in out_root else out_root

    # Determine if project is 2D or 3D
    is_2d_model, auto_halo = inspect_ilp(args.project)
    halo = args.halo if args.halo is not None else auto_halo

    n_classes = len(extract_ilp_labels(args.project))
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

    tiles = [{
            "y_min": args.y_min, "y_max": args.y_max,
            "x_min": args.x_min, "x_max": args.x_max,
            "z_min": args.z_min, "z_max": args.z_max,
            "t_index": t
            }
        for t in range(t_start, t_end)
        ]

    for tile in tiles:
        process_single_tile(tile, in_store, out_store, shell, args, is_2d_model, halo, n_classes)

    # Write single done file for the batch
    if args.done_file:
        done_path = Path(args.done_file)
        done_path.parent.mkdir(parents=True, exist_ok=True)
        with open(done_path, "w") as f:
            json.dump({"processed_count": len(tiles), "status": "COMPLETED", "tiles": tiles}, f, indent=2)


if __name__ == "__main__":
    main()
