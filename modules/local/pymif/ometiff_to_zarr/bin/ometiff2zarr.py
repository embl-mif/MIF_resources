#!/usr/bin/env python3

import argparse
import math
import sys
import pymif.microscope_manager as mm

def get_chunk_size(img_shape, chunk_size_input):
    """
    Calculates optimal chunk sizes to stay under 100MB per chunk.
    img_shape: [T, C, Z, Y, X]
    """
    MAX_BYTES = 100 * 1024 * 1024  # 100 MB
    BYTES_PER_PIXEL = 2            # 16-bit data

    if chunk_size_input == -1:
        # Prioritize spatial XY chunks for better visualization performance
        target_z = max(1, img_shape[2] // 4)
        target_y = min(img_shape[3], 4096)
        target_x = min(img_shape[4], 4096)

        size_bytes = target_z * target_y * target_x * BYTES_PER_PIXEL

        while size_bytes > MAX_BYTES:
            if target_x > 256 or target_y > 256:
                target_x = max(1, target_x // 2)
                target_y = max(1, target_y // 2)
            else:
                target_z = max(1, target_z // 2)
            size_bytes = target_z * target_y * target_x * BYTES_PER_PIXEL

        chunk_size_output = (1, 1, int(target_z), int(target_y), int(target_x))
    else:
        chunk_size_output = tuple(int(x) for x in chunk_size_input)

    n_chunks = [
        math.ceil(img_shape[2] / chunk_size_output[2]), 
        math.ceil(img_shape[3] / chunk_size_output[3]), 
        math.ceil(img_shape[4] / chunk_size_output[4])  
    ]

    return chunk_size_output, n_chunks

def get_num_layers(img_shape, num_layers_input):
    """Calculates pyramid levels."""
    z, y, x = img_shape[2], img_shape[3], img_shape[4]
    
    if num_layers_input == -1:
        num_layers_output = 1
        curr_z, curr_y, curr_x = z, y, x
        while curr_z > 2048 or curr_y > 2048 or curr_x > 2048:
            num_layers_output += 1
            curr_z //= 2
            curr_y //= 2
            curr_x //= 2
    else:
        num_layers_output = int(num_layers_input)

    return num_layers_output

def ometiff2zarr(input_ometiff, output_zarr, num_layers, chunk_size):
    # Initialize manager to read metadata
    dataset = mm.OperaManager(input_ometiff)
    img_shape = dataset.metadata["size"][0] # [T, C, Z, Y, X]

    # Calculate parameters
    final_chunks, _ = get_chunk_size(img_shape, chunk_size)
    final_layers = get_num_layers(img_shape, num_layers)

    print(f"Processing: {input_ometiff}")
    print(f"Metadata Shape: {img_shape}")
    print(f"Applied Chunks: {final_chunks} | Layers: {final_layers}")

    # Re-initialize with chunking and build
    dataset = mm.OperaManager(input_ometiff, chunks=final_chunks)
    dataset.build_pyramid(num_levels=final_layers, downscale_factor=2)
    dataset.to_zarr(output_zarr)
    print(f"Successfully created: {output_zarr}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert OME-TIFF to OME-ZARR using pymif.")
    parser.add_argument("-i", "--input_ometiff", required=True)
    parser.add_argument("-o", "--output_zarr", required=True)
    parser.add_argument("-nl", "--num_layers", type=str, default="")
    parser.add_argument("-cs", "--chunk_size", type=str, default="")

    args = parser.parse_args()

    # Parse nl
    nl = int(args.num_layers) if args.num_layers.strip() else -1

    # Parse cs (handles space-separated or comma-separated strings from Nextflow)
    if args.chunk_size.strip():
        raw_cs = args.chunk_size.replace(',', ' ').split()
        cs = tuple(int(x) for x in raw_cs)
    else:
        cs = -1

    ometiff2zarr(args.input_ometiff, args.output_zarr, nl, cs)
