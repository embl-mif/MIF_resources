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
    Parses OME-Zarr metadata, flexibly matching the intensity component.
    """
    reader = Reader(parse_url(input_zarr))
    nodes = list(reader())
    
    # Debug print to confirm exactly what we are seeing
    # print(f"DEBUG: Looking for '{component_path}' in nodes with paths: {[getattr(n.zarr, 'path', 'root') for n in nodes]}")

    for node in nodes:
        node_path = getattr(node.zarr, 'path', '')
        
        # 1. Exact match (standard case)
        if component_path == node_path:
            return node.data[0]
            
        # 2. Relaxed match: If we are looking for '0' and it's the root node
        # (detected by checking if it's the only node or contains the filename)
        if component_path == "0" and (node_path == "" or input_zarr in node_path):
            return node.data[0]
            
        # 3. Handle 'labels/nuclei' (path-based match)
        if component_path in node_path:
            return node.data[0]
            
    raise ValueError(f"Could not find component '{component_path}'. Available paths: {[getattr(n.zarr, 'path', 'root') for n in nodes]}")
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
        intensity_stack_all_ch = get_ome_zarr_array(input_zarr, intensity_path)
        intensity_stack = intensity_stack_all_ch[:, 0, :, :, :]# make arg later, defaults to 0 now
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
