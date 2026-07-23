#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import zarr
import vigra
from collections import OrderedDict

from ilastik import app
from ilastik.applets.dataSelection.opDataSelection import PreloadedArrayDatasetInfo


def _parse_optional_int(val):
    """Safely parse integer arguments that may be passed as empty strings from shell/Nextflow."""
    if val is None or val == "" or str(val).strip() == "":
        return None
    return int(val)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Ilastik pixel classification in-memory on a single 2D or 3D Zarr subregion with halo padding."
    )
    parser.add_argument("--input-zarr", required=True, help="Path to input raw Zarr store")
    parser.add_argument("--output-zarr", required=True, help="Path to preallocated output Zarr store")
    parser.add_argument("--project", required=True, help="Path to Ilastik project file (.ilp)")
    
    parser.add_argument("--t-index", type=int, default=0, help="Time index to process")

    # Z bounds are optional to seamlessly support 2D workflows
    # Flexible Z arguments: supports either a single slice (--z-index) or a range (--z-min/--z-max)
    parser.add_argument("--z-index", type=_parse_optional_int, default=None, help="Single Z slice index (for 2D model/mode)")
    parser.add_argument("--z-min", type=_parse_optional_int, default=None, help="Minimum Z index (optional for 2D)")
    parser.add_argument("--z-max", type=_parse_optional_int, default=None, help="Maximum Z index (optional for 2D)")
    
    parser.add_argument("--y-min", type=int, required=True, help="Minimum Y index")
    parser.add_argument("--y-max", type=int, required=True, help="Maximum Y index")
    parser.add_argument("--x-min", type=int, required=True, help="Minimum X index")
    parser.add_argument("--x-max", type=int, required=True, help="Maximum X index")
    parser.add_argument("--halo", type=int, default=16, help="Halo size in pixels")
    return parser.parse_args()

def check_ilp_is_2d(project_path):
    """Inspects the .ilp file to see if features are configured to compute in 2D."""
    with h5py.File(project_path, "r") as f:
        if "FeatureSelections/ComputeIn2d" in f:
            val = f["FeatureSelections/ComputeIn2d"][()]
            # If all or any are True, it's operating in 2D mode
            return np.all(val) or np.any(val)
    return True # Default fallback

def open_zarr_array(path, mode="r"):
    """Open Zarr store, automatically selecting scale '0' if it is an OME-Zarr group. Only works on Zarr v2, OME-Zarr 0.4"""
    store = zarr.open(path, mode=mode)
    if isinstance(store, zarr.hierarchy.Group) and "0" in store:
        return store["0"]
    return store

def main():
    args = parse_args()
    # Resource management for Ilastik's lazyflow engine
    os.environ["LAZYFLOW_THREADS"] = os.environ.get("LAZYFLOW_THREADS", "2")
    os.environ["LAZYFLOW_TOTAL_RAM_MB"] = os.environ.get("LAZYFLOW_TOTAL_RAM_MB", "4000")

    # Open Zarr stores
    in_store = open_zarr_array(args.input_zarr, mode="r")
    out_store = open_zarr_array(args.output_zarr, mode="r+")

    # Determine if project is 2D or 3D
    is_2d_model = check_ilp_is_2d(args.project)
    
    # Determine Z operational mode
    # If explicit z-min/z-max are passed, treat as a 3D block job
    is_3d_block = (args.z_min is not None) and (args.z_max is not None)

    shape_t, shape_c, shape_z, shape_y, shape_x = in_store.shape

    # --- 3D MODE (Volumetric block with Z halo) ---
        pz_min = max(0, args.z_min - args.halo)
        pz_max = min(shape_z, args.z_max + args.halo)
        z_read_slice = slice(pz_min, pz_max)
        
        cz_start = args.halo if args.z_min > 0 else 0
        cz_end = cz_start + (args.z_max - args.z_min)
        z_crop_slice = slice(cz_start, cz_end)
        z_write_slice = slice(args.z_min, args.z_max)
    else:
        # --- 2D MODE (Single slice, no Z halo needed) ---
        z_idx = args.z_index if args.z_index is not None else 0
        z_read_slice = slice(z_idx, z_idx + 1)
        z_crop_slice = slice(0, 1)
        z_write_slice = slice(z_idx, z_idx + 1)

    py_min = max(0, args.y_min - args.halo)
    py_max = min(shape_y, args.y_max + args.halo)
    px_min = max(0, args.x_min - args.halo)
    px_max = min(shape_x, args.x_max + args.halo)

    spatial_read_slices.extend([slice(py_min, py_max), slice(px_min, px_max)])

    # Construct strict 5D read slice for [t, c, z, y, x]
    in_read_slices = (
        slice(args.t_index, args.t_index + 1),  # t
        slice(0, 1),                            # c
        z_read_slice,                           # z
        slice(py_min, py_max),                  # y
        slice(px_min, px_max)                   # x
    )

    # Read raw chunk + halo into RAM
    raw_chunk = in_store[tuple(in_read_slices)]

    # Initialize Ilastik
    ilastik_args = app.parse_args([])
    ilastik_args.headless = True
    ilastik_args.project = args.project
    
    shell = app.main(ilastik_args)
    
    # Tag axes for Ilastik
    axistags = "tczyx"
    tagged_chunk = vigra.taggedView(raw_chunk, axistags)
    role_data_dict = OrderedDict([
        ("Raw Data", [PreloadedArrayDatasetInfo(preloaded_array=tagged_chunk)]),
    ])

    # Execute prediction
    predictions = shell.workflow.batchProcessingApplet.run_export(role_data_dict, export_to_array=True)
    prob_chunk = predictions[0]

    # Align channel axis back to position 1 for TCZYX layout
    if prob_chunk.ndim == 5 and prob_chunk.shape[-1] == out_store.shape[1]:
        prob_chunk = np.moveaxis(prob_chunk, -1, 1)

    # Calculate halo cropping offsets relative to chunk shape
    cy_start = args.halo if args.y_min > 0 else 0
    cy_end = cy_start + (args.y_max - args.y_min)
    cx_start = args.halo if args.x_min > 0 else 0
    cx_end = cx_start + (args.x_max - args.x_min)

    # Crop halo from prediction output [t, c, z, y, x]
    crop_slices = (
        slice(None),  # t
        slice(None),  # c
        z_crop_slice,
        slice(cy_start, cy_end),
        slice(cx_start, cx_end)
    )
    cropped_probs = prob_chunk[crop_slices]

    # Convert probabilities (0.0 to 1.0) to uint8 (0 to 255) if destination store is uint8
    if out_store.dtype == np.uint8 and cropped_probs.dtype != np.uint8:
        cropped_probs = (np.clip(cropped_probs, 0, 1) * 255).astype(np.uint8)

    # Write directly into the allocated Zarr store
    write_slices = (
        slice(args.t_index, args.t_index + 1),
        slice(None),  # All probability channels
        z_write_slice,
        slice(args.y_min, args.y_max),
        slice(args.x_min, args.x_max)
    )
    out_store[write_slices] = cropped_probs
    print(f"Successfully populated T:{args.t_index}, Z-mode:{'3D-block' if is_3d_block else '2D-slice'}, Y:{args.y_min}-{args.y_max}, X:{args.x_min}-{args.x_max}")


if __name__ == "__main__":
    main()
