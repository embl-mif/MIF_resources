#!/usr/bin/env python3

import os
import glob
import argparse
import csv
import pandas as pd
import tifffile
import xml.etree.ElementTree as ET

def get_ome_metadata(ome_tif_path):
    """Extracts microscope model, channel names, and multiscale levels from an OME-TIFF."""
    meta = {
        "microscope": "unknown",
        "channel_names": "",
        "num_levels": ""
    }
    
    try:
        with tifffile.TiffFile(ome_tif_path) as tif:
            # 1. Extract pyramid levels directly from the TIFF hierarchy
            if tif.series:
                # series[0].levels contains the base image + all downsampled sub-resolutions
                meta["num_levels"] = str(len(tif.series[0].levels))
            
            # 2. Parse the OME-XML string
            ome_xml = tif.ome_metadata
            if not ome_xml: 
                return meta
                
            root = ET.fromstring(ome_xml)
            
            # 3. Extract and map the Microscope model
            # Using {*} wildcard to avoid strict OME namespace version dependencies
            microscope = root.find('.//{*}Microscope')
            if microscope is not None:
                model = microscope.get('Model', '').lower()
                # Force mapping to 'opera' for Phenix systems to satisfy pymif
                if 'phenix' in model or 'opera' in model:
                    meta["microscope"] = "opera"
                else:
                    meta["microscope"] = model

            # 4. Extract Channel Names
            # Look at the first image's pixels to find the channels
            first_pixels = root.find('.//{*}Pixels')
            if first_pixels is not None:
                channels = first_pixels.findall('./{*}Channel')
                ch_names = []
                for ch in channels:
                    name = ch.get('Name')
                    if name:
                        ch_names.append(name)
                    else:
                        # Fallback if 'Name' is missing
                        ch_names.append(ch.get('ID', 'unknown'))
                
                # Join with commas (e.g., "DAPI,Alexa 488,Hoechst")
                meta["channel_names"] = ",".join(ch_names)

    except Exception as e:
        print(f"Warning: Failed to parse metadata for {ome_tif_path}: {e}")
        
    return meta

def generate_sample_input(exp_folder, zarr_version, output_csv):
    # Construct the full path to the data directory
    
    if not os.path.exists(exp_folder):
        raise RuntimeError(f"Directory {exp_folder} does not exist.")

    headers = ["input", "microscope", "output", "chunk_size", "max_size(MB)", 
               "scene_index", "zarr_format", "downscale_factor", "channel_colors", 
               "channel_names", "num_levels"]

    meta_paths = glob.glob(f"{exp_folder}/*metadata.csv")
    samples_to_process = []

    if meta_paths:
        print(f"Metadata found: {meta_paths[0]}. Filtering samples...")
        df = pd.read_csv(meta_paths[0])
        if "ignore" in df.columns:
            df = df[~df["ignore"].astype(str).str.lower().isin(["true", "1", "yes"])]
        samples_to_process = df["sample"].tolist()
    else:
        print("No metadata found. Processing all OME-TIFFs in directory...")
        all_files = glob.glob(f"{exp_folder}/*.ome.tiff") + glob.glob(f"{exp_folder}/*.ome.tif")
        samples_to_process = [os.path.basename(f).split('.ome')[0] for f in all_files]

    processed_count = 0

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for sample in samples_to_process:
            # Find files starting with the sample name in the target_dir
            ome_files = glob.glob(f"{exp_folder}/{sample}*.ome.tiff") + glob.glob(f"{exp_folder}/{sample}*.ome.tif")
            if not ome_files: continue
            
            ome_file = ome_files[0]
            # Extract metadata for this specific file
            file_meta = get_ome_metadata(ome_file)
            
            writer.writerow([
                ome_file, file_meta["microscope"], f"{sample}.zarr", "1 1 1 512 512",
                "100", "0", args.zarr_version, "2", "", file_meta["channel_names"], file_meta["num_levels"]
            ])
            processed_count += 1

    print(f"Sample input table complete. Wrote {processed_count} entries to {output_csv}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", required=True, help="Path to the experiment root.")
    parser.add_argument("--zarr-version", default="2", help="Zarr version to use (2 or 3). Default: 2")
    parser.add_argument("--out", default="sample_input.csv")
    args = parser.parse_args()
    generate_sample_input(args.exp, args.zarr_version, args.out)
