#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
import dask
import dask.array as da
import pandas as pd
from skimage.measure import regionprops_table
import pyarrow

# Configure logging for Nextflow compatibility
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)

@dask.delayed
def process_plane(intensity_plane, label_plane, z_index):
    """
    Processes a single Z-slice using skimage regionprops.
    Assumes intensity_plane is (C, Y, X).
    """
    # Convert CZYX -> YXC for skimage
    # Adjust this transposition if your intensity data has a different axis order
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
        logging.info(f"Opening Zarr store: {input_zarr}")
        
        # Load Zarr groups lazily as Dask arrays
        intensity_stack = da.from_zarr(input_zarr, component=intensity_path)
        label_stack = da.from_zarr(input_zarr, component=label_path)

        if intensity_stack.shape[1:3] != label_stack.shape[1:3]:
            raise ValueError(f"Shape mismatch: Intensity {intensity_stack.shape} vs Labels {label_stack.shape}")

        n_planes = label_stack.shape[0]
        logging.info(f"Found {n_planes} planes. Launching Dask tasks...")

        # Create delayed tasks
        tasks = []
        for i in range(n_planes):
            tasks.append(process_plane(
                intensity_stack[:, i, ...], 
                label_stack[i, ...], 
                i
            ))

        logging.info("Computing parallel measurements...")
        results = dask.compute(*tasks) 

        # Combine results
        final_df = pd.concat(results, ignore_index=True)
        
        # Save as Parquet
        final_df.to_parquet(output_path, index=False, engine='pyarrow')
        logging.info(f"Successfully saved output to {output_path}")

    except Exception as e:
        logging.error(f"Error processing {input_zarr}: {str(e)}")
        sys.exit(1) # Ensure Nextflow receives non-zero exit code on failure

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from Zarr for Nextflow.")
    parser.add_argument("-z", "--input_zarr", type=str, required=True, help="Input Zarr path")
    parser.add_argument("-o", "--output_path", type=str, required=True, help="Output Parquet path")
    parser.add_argument("--intensity", type=str, default="0", help="Zarr path for intensity data")
    parser.add_argument("--labels", type=str, default="labels/nuclei", help="Zarr path for labels")

    args = parser.parse_args()

    run_measurement(args.input_zarr, args.output_path, args.intensity, args.labels)
