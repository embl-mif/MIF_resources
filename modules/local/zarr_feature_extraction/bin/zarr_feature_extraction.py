#!/usr/bin/env python3
import argparse
import logging
import sys
import dask
import dask.array as da
import pandas as pd
from skimage.measure import regionprops_table
import pyarrow

# OME-Zarr imports for V5 compatibility
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)

def get_ome_zarr_array(input_zarr, component_path):
    """
    Parses OME-Zarr V5 metadata to find the specific array node.
    """
    reader = Reader(parse_url(input_zarr))
    # Iterate through nodes to find the one matching your path/component
    for node in reader():
        # component_path usually looks like '0' or 'labels/nuclei'
        if component_path in node.path:
            # node.data is a list of dask arrays (pyramid levels)
            # node.data[0] is typically the full-resolution image
            return node.data[0]
    
    raise ValueError(f"Could not find component '{component_path}' in OME-Zarr V5 store.")

@dask.delayed
def process_plane(intensity_plane, label_plane, z_index):
    # Convert CZYX -> YXC for skimage
    intensity_img = intensity_plane.transpose(1, 2, 0)
    
    table = regionprops_table(
        label_image=label_plane, 
        intensity_image=intensity_img, 
        properties=[
            "label", "area", "centroid", "intensity_mean", 
            "intensity_max", "intensity_min", "intensity_std", "moments"
        ]
    )

    df = pd.DataFrame(table)
    df["z"] = z_index
    return df

def run_measurement(input_zarr, output_path, intensity_path, label_path):
    try:
        logging.info(f"Opening OME-Zarr V5 store: {input_zarr}")
        
        # Use the specialized reader instead of da.from_zarr
        intensity_stack = get_ome_zarr_array(input_zarr, intensity_path)
        label_stack = get_ome_zarr_array(input_zarr, label_path)

        # Ensure shapes match (using Y, X, C order)
        if intensity_stack.shape[1:3] != label_stack.shape[1:3]:
            raise ValueError(f"Shape mismatch: Intensity {intensity_stack.shape} vs Labels {label_stack.shape}")

        n_planes = label_stack.shape[0]
        logging.info(f"Found {n_planes} planes. Launching Dask tasks...")

        tasks = []
        for i in range(n_planes):
            tasks.append(process_plane(
                intensity_stack[:, i, ...], 
                label_stack[i, ...], 
                i
            ))

        logging.info("Computing parallel measurements...")
        results = dask.compute(*tasks) 

        final_df = pd.concat(results, ignore_index=True)
        final_df.to_parquet(output_path, index=False, engine='pyarrow')
        logging.info(f"Successfully saved output to {output_path}")

    except Exception as e:
        logging.error(f"Error processing {input_zarr}: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from OME-Zarr V5 for Nextflow.")
    parser.add_argument("-z", "--input_zarr", type=str, required=True, help="Input Zarr path")
    parser.add_argument("-o", "--output_path", type=str, required=True, help="Output Parquet path")
    parser.add_argument("--intensity", type=str, default="0", help="Component for intensity (e.g., '0')")
    parser.add_argument("--labels", type=str, default="labels/nuclei", help="Component for labels (e.g., 'labels/nuclei')")

    args = parser.parse_args()
    run_measurement(args.input_zarr, args.output_path, args.intensity, args.labels)
