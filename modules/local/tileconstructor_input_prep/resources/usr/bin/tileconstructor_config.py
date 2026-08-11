#!/usr/bin/env python3
import os
import json
import csv
import argparse

def get_zarr_metadata(zarr_path, target_level):
    zarr_path = zarr_path.rstrip('/')
    is_v3 = os.path.exists(os.path.join(zarr_path, "zarr.json"))
    
    if is_v3:
        with open(os.path.join(zarr_path, "zarr.json"), "r") as f:
            root_meta = json.load(f)
            
        datasets = root_meta.get("multiscales", [{}])[0].get("datasets", [])
        level_path = next((ds.get("path") for ds in datasets if str(ds.get("path")) == str(target_level)), None)
        if level_path is None:
            raise ValueError(f"Resolution level {target_level} not found in {zarr_path}")
            
        level_json_path = os.path.join(zarr_path, str(level_path), "zarr.json")
        with open(level_json_path, "r") as f:
            level_meta = json.load(f)
            
        shape = level_meta.get("shape")
        chunk_config = level_meta.get("chunk_grid", {}).get("configuration", {}).get("chunk_shape", [])
        chunks = chunk_config if chunk_config else shape
        return shape, chunks
    else:
        level_array_path = os.path.join(zarr_path, str(target_level), ".zarray")
        if not os.path.exists(level_array_path):
            raise ValueError(f"Metadata file not found at {level_array_path}")
            
        with open(level_array_path, "r") as f:
            level_meta = json.load(f)
            
        shape = level_meta.get("shape")
        chunks = level_meta.get("chunks")
        return shape, chunks

def scan_directory_for_ome_zarr(root_directory, target_level=0):
    datasets_info = []
    if not os.path.exists(root_directory):
        raise NotADirectoryError(f"Directory not found: {root_directory}")

    root_directory = os.path.abspath(root_directory)
    for entry in os.listdir(root_directory):
        full_path = os.path.join(root_directory, entry)
        if os.path.isdir(full_path) and entry.endswith(".zarr"):
            is_v3 = os.path.exists(os.path.join(full_path, "zarr.json"))
            is_v2 = os.path.exists(os.path.join(full_path, ".zgroup")) or os.path.exists(os.path.join(full_path, ".zarray"))
            if is_v3 or is_v2:
                datasets_info.append({
                    "path": full_path,
                    "level": target_level,
                    "dataset_id": f"{os.path.splitext(entry)[0]}"
                })
    return datasets_info

def generate_tile_config_csv(root_directory, output_csv="tile_constructor_config.csv", target_level=0):
    datasets_info = scan_directory_for_ome_zarr(root_directory, target_level)
    headers = [
        "dataset_id", "input_image_path", "tile_size_x", "tile_size_y", 
        "tile_size_z", "tile_overlap", "resolution_level", 
        "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"
    ]
    
    rows = []
    for info in datasets_info:
        zarr_path = info["path"]
        level = info["level"]
        dataset_id = info["dataset_id"]
        
        shape, chunks = get_zarr_metadata(zarr_path, level)
        
        if len(shape) >= 5 and len(chunks) >= 5:
            z_max = shape[-3]
            y_max = shape[-2]
            x_max = shape[-1]
            tile_z = chunks[-3]
            tile_y = chunks[-2]
            tile_x = chunks[-1]
        else:
            raise ValueError(f"Unexpected shape/chunk format in {zarr_path}")
            
        overlap = 0  
        
        rows.append([
            dataset_id, zarr_path, tile_x, tile_y, tile_z, overlap, level,
            0, x_max, 0, y_max, 0, z_max
        ])
        
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
        
    print(f"Successfully generated configuration CSV: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-csv", default="tile_constructor_config.csv")
    parser.add_argument("--resolution-level", type=int, default=0)
    args = parser.parse_args()
    
    generate_tile_config_csv(args.input_dir, args.output_csv, args.resolution_level)
