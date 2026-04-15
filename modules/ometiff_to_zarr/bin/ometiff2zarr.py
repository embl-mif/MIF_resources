#!/usr/bin/env python3

import argparse
import math
import pymif.microscope_manager as mm

def get_chunk_size(img_shape, chunk_size_input):
    MAX_BYTES = 100 * 1024 * 1024  # 100 MB
    BYTES_PER_PIXEL = 2            # 16-bit data

    if chunk_size_input == -1:
        # Start with full image for X and Y, and 1 for Z
        # We prioritize keeping X and Y large for spatial performance
        target_z = max(1, img_shape[2] // 4)
        target_y = min(img_shape[3], 4096)
        target_x = min(img_shape[4], 4096)
        
        # Calculate current size
        size_bytes = target_z * target_y * target_x * BYTES_PER_PIXEL
            
        # If still too big, keep shrinking X and Y by half
        while size_bytes > MAX_BYTES:
            if target_x > 256 or target_y > 256:
                target_x = max(1, target_x // 2)
                target_y = max(1, target_y // 2)
            else:
                # If xy chunks is already small, shrink Z
                target_z = max(1, target_z // 2)
            size_bytes = target_z * target_y * target_x * BYTES_PER_PIXEL

            if size_bytes > MAX_BYTES:
                target_z = max(1, target_z // 2)
                size_bytes = target_z * target_y * target_x * BYTES_PER_PIXEL
            
        chunk_size_output = (1, 1, int(target_z), int(target_y), int(target_x))
        
    else:
        # 1. Ensure it's a tuple of ints
        chunk_size_output = tuple(int(x) for x in chunk_size_input)

    size_bytes = chunk_size_output[2] * chunk_size_output[3] * chunk_size_output[4] * BYTES_PER_PIXEL

    # 2. Compute n_chunks based on the user-provided chunk_size
    # Formula: ceil(Total_Size / Chunk_Size)
    n_chunks = [
        math.ceil(img_shape[2] / chunk_size_output[2]), # Z
        math.ceil(img_shape[3] / chunk_size_output[3]), # Y
        math.ceil(img_shape[4] / chunk_size_output[4])  # X
    ]

    return chunk_size_output, n_chunks, size_bytes

def get_num_layers(img_shape, num_layers_input):

    shape = [[img_shape[2], img_shape[3], img_shape[4]]] # [Z, Y, X]

    if num_layers_input == -1:

        num_layers_output = 1
        # print(f"Layer {n}, shape {shape}")

        while (shape[-1][0]>2048) or (shape[-1][1]>2048) or (shape[-1][2]>2048):

            num_layers_output+=1
            # print(f"Layer {n}, shape {shape}")
            shape.append( [shape[-1][0]//2, shape[-1][1]//2, shape[-1][2]//2] )
    else:
        num_layers_output = int(num_layers_input)

        for n in range(1, num_layers_output):
            shape.append( [shape[-1][0]//2, shape[-1][1]//2, shape[-1][2]//2] )
            # print(f"Layer {n}, shape {shape}")

    return num_layers_output, shape

def ometiff2zarr(input_ometiff, output_zarr, num_layers=3, chunk_size=(1,1,8,256,256)):

    dataset = mm.OperaManager(input_ometiff)

    print(dataset.metadata)

    # 1. Handle Chunk Size Logic
    # Assume dataset.metadata["size"][0] is [T, C, Z, Y, X]
    img_shape = dataset.metadata["size"][0]

    chunk_size, n_chunks, size_bytes = get_chunk_size(img_shape, chunk_size)

    print(f"Final Chunk Size: {chunk_size}")
    print(f"Grid Layout (Z,Y,X chunks): {n_chunks}")
    print(f"Size in Bytes: {size_bytes}")

    # 2. Handle Num Layers Logic

    num_layers, shapes = get_num_layers(img_shape, num_layers)

    print(f"Final Num Layers: {num_layers}")
    print(f"Layer shapes: {shapes}")
    
    print(f"Final params -> Chunks: {chunk_size}, Layers: {num_layers}")
    
    dataset = mm.OperaManager(input_ometiff, chunks=chunk_size)
    dataset.build_pyramid(num_levels=num_layers, downscale_factor=2)
    dataset.to_zarr(output_zarr)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_ometiff", type=str, required=True)
    parser.add_argument("-o", "--output_zarr", type=str, required=True)
    
    # Use strings for optional inputs so we can detect "empty" from Nextflow
    parser.add_argument('-nl', '--num_layers', type=str, default="")
    parser.add_argument('-cs', '--chunk_size', type=str, default="")

    args = parser.parse_args()

    # Determine num_layers
    if not args.num_layers or args.num_layers == "":
        nl = -1
    else:
        nl = int(args.num_layers)

    # Determine chunk_size
    if not args.chunk_size or args.chunk_size == "":
        cs = -1
    else:
        # Pass the string to the function; the function handles the split
        cs = tuple(int(x.strip()) for x in str(args.chunk_size).split(' '))

    ometiff2zarr(args.input_ometiff, args.output_zarr,
                 num_layers=nl, chunk_size=cs)
