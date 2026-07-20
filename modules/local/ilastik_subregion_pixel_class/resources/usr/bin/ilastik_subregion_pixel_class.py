#!/usr/bin/env python3
import argparse
import sys
import numpy as np
import zarr
from ilastik.applets.pixelClassification import PixelClassificationHeadlessApplet
from ilastik.shell import IlastikShell


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
    
    # Z bounds are optional to seamlessly support 2D workflows
    parser.add_argument("--z-min", type=_parse_optional_int, default=None, help="Minimum Z index (optional for 2D)")
    parser.add_argument("--z-max", type=_parse_optional_int, default=None, help="Maximum Z index (optional for 2D)")
    
    parser.add_argument("--y-min", type=int, required=True, help="Minimum Y index")
    parser.add_argument("--y-max", type=int, required=True, help="Maximum Y index")
    parser.add_argument("--x-min", type=int, required=True, help="Minimum X index")
    parser.add_argument("--x-max", type=int, required=True, help="Maximum X index")
    parser.add_argument("--halo", type=int, default=16, help="Halo size in pixels")
    return parser.parse_args()


def open_zarr_array(path, mode="r"):
    """Open Zarr store, automatically selecting scale '0' if it is an OME-Zarr group. Only works on Zarr v2, OME-Zarr 0.4"""
    store = zarr.open(path, mode=mode)
    if isinstance(store, zarr.hierarchy.Group) and "0" in store:
        return store["0"]
    return store


def main():
    args = parse_args()

    # 1. Open Zarr stores
    in_store = open_zarr_array(args.input_zarr, mode="r")
    out_store = open_zarr_array(args.output_zarr, mode="r+")

    # Determine if this task is operating in 3D or 2D mode
    has_z_args = (args.z_min is not None) and (args.z_max is not None)
    
    # Check if the store actually has enough spatial dimensions for Z
    # (Assuming spatial dimensions occupy trailing axes)
    if has_z_args and in_store.ndim >= 3 and in_store.shape[-3] > 1:
        has_z = True
    else:
        has_z = False
        if has_z_args and in_store.ndim < 3:
            print("Note: --z-min/--z-max provided, but input data is 2D. Running in 2D mode.")

    # 2. Extract spatial shapes and calculate padded bounds with halo
    if has_z:
        shape_z, shape_y, shape_x = in_store.shape[-3:]
        pz_min = max(0, args.z_min - args.halo)
        pz_max = min(shape_z, args.z_max + args.halo)
        spatial_read_slices = [slice(pz_min, pz_max), slice(None), slice(None)]
    else:
        shape_y, shape_x = in_store.shape[-2:]
        spatial_read_slices = []

    py_min = max(0, args.y_min - args.halo)
    py_max = min(shape_y, args.y_max + args.halo)
    px_min = max(0, args.x_min - args.halo)
    px_max = min(shape_x, args.x_max + args.halo)

    spatial_read_slices.extend([slice(py_min, py_max), slice(px_min, px_max)])

    # Construct dynamic slice object for reading (filling non-spatial leading dims with :)
    num_spatial_dims = 3 if has_z else 2
    in_read_slices = [slice(None)] * (in_store.ndim - num_spatial_dims) + spatial_read_slices

    # 3. Read raw chunk + halo into RAM
    raw_chunk = in_store[tuple(in_read_slices)]

    # 4. Initialize Ilastik Headless Applet in RAM
    shell = IlastikShell()
    applet = PixelClassificationHeadlessApplet()
    op_pixel = applet.topLevelOperator
    op_pixel.ProjectFile.setValue(args.project)

    # Set input data and execute prediction pipeline
    op_pixel.InputImages.setValues([raw_chunk])
    prob_chunk = op_pixel.HeadlessPredictionPipeline[0].export_data().wait()

    # 5. Calculate halo cropping offsets
    if has_z:
        cz_start = args.halo if args.z_min > 0 else 0
        cz_end = cz_start + (args.z_max - args.z_min)
        spatial_crop_slices = [slice(cz_start, cz_end)]
        spatial_write_slices = [slice(args.z_min, args.z_max)]
    else:
        spatial_crop_slices = []
        spatial_write_slices = []

    cy_start = args.halo if args.y_min > 0 else 0
    cy_end = cy_start + (args.y_max - args.y_min)
    cx_start = args.halo if args.x_min > 0 else 0
    cx_end = cx_start + (args.x_max - args.x_min)

    spatial_crop_slices.extend([slice(cy_start, cy_end), slice(cx_start, cx_end)])
    spatial_write_slices.extend([slice(args.y_min, args.y_max), slice(args.x_min, args.x_max)])

    # Crop halo from prediction output
    crop_slices = [slice(None)] * (prob_chunk.ndim - len(spatial_crop_slices)) + spatial_crop_slices
    cropped_probs = prob_chunk[tuple(crop_slices)]

    # Convert probabilities (0.0 to 1.0) to uint8 (0 to 255) if destination store is uint8
    if out_store.dtype == np.uint8 and cropped_probs.dtype != np.uint8:
        cropped_probs = (np.clip(cropped_probs, 0, 1) * 255).astype(np.uint8)

    # 6. Write directly into the allocated Zarr store
    write_slices = [slice(None)] * (out_store.ndim - len(spatial_write_slices)) + spatial_write_slices
    out_store[tuple(write_slices)] = cropped_probs

    z_str = f"Z:{args.z_min}-{args.z_max}, " if has_z else ""
    print(f"Successfully populated slice {z_str}Y:{args.y_min}-{args.y_max}, X:{args.x_min}-{args.x_max}")


if __name__ == "__main__":
    main()
