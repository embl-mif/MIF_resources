#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
import dask
from dask import delayed
import pandas as pd
import pymif.microscope_manager as mm
from skimage.measure import regionprops_table

# Configure logging for Nextflow compatibility
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stdout
)

@delayed
def process_plane(intensity_plane, label_plane, z_index):
    """Processes a single Z-slice."""
    # Convert CZYX -> YXC for skimage
    intensity_img = intensity_plane.transpose(1, 2, 0)
    table = regionprops_table(
        label_image=label_plane, 
        intensity_image=intensity_img, 
        properties=["label", "area", "centroid", "intensity_mean", 
                    "intensity_max", "intensity_min", "intensity_std", "moments"]
    )

    df = pd.DataFrame(table)
    df["z"] = z_index
    return df

def run_measurement(input_zarr, output_csv):
    try:
        logging.info(f"Loading dataset: {input_zarr}")
        dataset = mm.ZarrManager(input_zarr)

        # Accessing the dask arrays
        intensity_stack = dataset.data[0] # Lazy
        label_stack = dataset.labels["nuclei"][0] # Lazy

        n_planes = label_stack.shape[0]
        logging.info(f"Found {n_planes} planes to process.")

        tasks = []
        for i in range(n_planes):
            # Pass dask slices, lazy
            tasks.append(process_plane(intensity_stack[:, i, ...], label_stack[i, ...], i))

        logging.info("Starting parallel feature extraction...")
        results = dask.compute(*tasks) # Load to memory here

        # Combine
        final_df = pd.concat(results, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        logging.info(f"Successfully saved output to {output_csv}")

    except Exception as e:
        logging.error(f"Failed to process {input_zarr}: {str(e)}")
        sys.exit(1) # Ensure Nextflow knows the process failed

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from Zarr for Nextflow.")
    parser.add_argument("-z", "--input_zarr", type=Path, required=True, help="Input Zarr path")
    parser.add_argument("-o", "--output_csv", type=Path, required=True, help="Output CSV path")

    args = parser.parse_args()

    run_measurement(str(args.input_zarr), str(args.output_csv))
