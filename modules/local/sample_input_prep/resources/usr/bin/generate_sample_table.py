#!/usr/bin/env python3

import os
import glob
import argparse
import csv
import pandas as pd
import tifffile
import xml.etree.ElementTree as ET

def get_microscope_name(ome_tif_path):
    try:
        with tifffile.TiffFile(ome_tif_path) as tif:
            ome_xml = tif.ome_metadata
            if not ome_xml: return "unknown"
            root = ET.fromstring(ome_xml)
            instruments = root.findall('.//{http://www.openmicroscopy.org/Schemas/OME/2016-06}Instrument')
            for instr in instruments:
                microscope = instr.find('./{http://www.openmicroscopy.org/Schemas/OME/2016-06}Microscope')
                if microscope is not None: return microscope.get('Model', 'unknown_model')
    except Exception: return "unknown"
    return "unknown"

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

    detected_microscope = None
    processed_count = 0

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for sample in samples_to_process:
            # Find files starting with the sample name in the target_dir
            ome_files = glob.glob(f"{exp_folder}/{sample}*.ome.tiff") + glob.glob(f"{exp_folder}/{sample}*.ome.tif")
            if not ome_files: continue
            
            ome_file = ome_files[0]
            if detected_microscope is None:
                detected_microscope = get_microscope_name(ome_file)
            
            writer.writerow([
                ome_file, detected_microscope, f"{sample}.zarr", "1 1 1 512 512",
                "100", "0", args.zarr_version, "2", "", "", ""
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
