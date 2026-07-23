#!/usr/bin/env python3
import argparse
import os
import shutil
import h5py
import zarr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preallocate an empty OME-Zarr skeleton for parallel writes. Use Ilastik .ilp project metadata"
    )
    parser.add_argument(
        "--input", required=True, help="Path to the input OME-Zarr directory"
    )
    parser.add_argument(
        "--ilp", required=True, help="Path to the Ilastik .ilp project file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path where the preallocated OME-Zarr will be created",
    )
    return parser.parse_args()

def extract_ilp_labels(ilp_path):
    """Extract label names and count directly from the Ilastik project file. Label # determines channel #"""
    if not os.path.exists(ilp_path):
        raise FileNotFoundError(f"Ilastik project file not found: {ilp_path}")

    with h5py.File(ilp_path, "r") as f:
        if "PixelClassification/LabelNames" not in f:
            raise KeyError(
                "Could not find '/PixelClassification/LabelNames' in the .ilp file."
            )
        
        raw_labels = f["PixelClassification/LabelNames"][()]
        labels = [
            l.decode("utf-8") if isinstance(l, bytes) else str(l)
            for l in raw_labels
        ]
        return label

def preallocate_zarr(input_path, ilp_path, output_path):
    # Read class metadata from Ilastik .ilp file
    labels = extract_ilp_labels(ilp_path)
    num_channels = len(labels)
    print(f"Detected {num_channels} classes from .ilp: {labels}")

    # Open the source Zarr group
    src_group = zarr.open(input_path, mode="r")

    # Create the destination Zarr directory
    if os.path.exists(output_path):
        print(f"Warning: Output path '{output_path}' already exists. Overwriting...")
        shutil.rmtree(output_path)
    
    # Use DirectoryStore explicitly for better parallel write compatibility
    store = zarr.DirectoryStore(output_path)
    dest_group = zarr.group(store=store, overwrite=True)

    # Read and modify OME-NGFF metadata (.zattrs)
    zattrs = dict(src_group.attrs)
    multiscales = zattrs.get("multiscales", [])
    if not multiscales:
        raise ValueError(f"Input Zarr '{input_path}' lacks 'multiscales' in .zattrs.")

    # Detect the coordinate system axes to find 'c' (channel)
    axes = multiscales[0].get("axes", [])
    channel_axis_index = None
    for idx, axis in enumerate(axes):
        axis_name = axis.get("name") if isinstance(axis, dict) else axis
        if axis_name == "c":
            channel_axis_index = idx
            break

    if channel_axis_index is None:
        channel_axis_index = 1
        print("Warning: Could not explicitly find channel axis 'c'. Defaulting to index 1.")

    # Update metadata to reflect the new channel count for viewers
    if "omero" in zattrs:
        # Standard distinct hex colors for microscopy/segmentation classes
        palette = ["FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "FFFFFF"]
        new_omero_channels = []
        
        for i, label_name in enumerate(labels):
            color = palette[i % len(palette)]
            new_omero_channels.append(
                {
                    "active": True,
                    "coefficient": 1.0,
                    "color": color,
                    "family": "linear",
                    "inverted": False,
                    "label": f"Probability Class {i}",
                    "window": {"end": 1.0, "max": 1.0, "min": 0.0, "start": 0.0},
                }
            )
        zattrs["omero"]["channels"] = new_omero_channels

    dest_group.attrs.update(zattrs)

    # Create each empty multiscale pyramid level array
    for dataset_meta in multiscales[0]["datasets"]:
        path_name = dataset_meta["path"]
        src_array = src_group[path_name]

        orig_shape = list(src_array.shape)
        orig_chunks = list(src_array.chunks)
        
        # Guard against malformed dimension mapping
        if channel_axis_index >= len(orig_shape):
            raise IndexError(
                f"Channel index {channel_axis_index} is out of bounds for shape {orig_shape}."
            )

        # Modify channel size
        new_shape = list(orig_shape)
        new_shape[channel_axis_index] = num_channels

        new_chunks = list(orig_chunks)
        new_chunks[channel_axis_index] = 1

        print(f"Preallocating scale '{path_name}': {orig_shape} -> {new_shape}")

        # Create empty array
        dest_group.create_dataset(
            name=path_name,
            shape=new_shape,
            chunks=new_chunks,
            dtype="float32", 
            compressor=src_array.compressor,
            fill_value=0.0,
        )

    print(f"Successfully preallocated empty OME-Zarr store at: {output_path}")


if __name__ == "__main__":
    args = parse_args()
    preallocate_zarr(
        input_path=args.input,
        ilp_path=args.ilp,
        output_path=args.output
    )
