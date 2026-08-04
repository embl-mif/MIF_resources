#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

from squirrel.library.ome_zarr import OMEZarrStore


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild OME-Zarr pyramid levels 1..N from level 0, using squirrel's PyramidBuilder."
    )
    parser.add_argument("--output-zarr", required=True,
                         help="Path to the OME-Zarr store to rebuild (level 0 must already be fully written).")
    parser.add_argument("--downsample-method", choices=["Average", "Sample"], default="Average",
                         help="Downsampling method. Use 'Average' for continuous data like probabilities, "
                              "'Sample' for label/instance images.")
    parser.add_argument("--n-threads", type=int, default=1, help="Thread count for pyramid rebuild.")

    args = parser.parse_args()

    print(f"[PyramidRebuild] Opening store: {args.output_zarr}")
    oz = OMEZarrStore(args.output_zarr, mode="a")

    # Set explicitly here because it's not parsing Zarr stores generated elsewere outside ome_zarr.py correctly
    oz.metadata.downsample_method = args.downsample_method

    print(f"[PyramidRebuild] levels: {oz.metadata.levels}")
    print(f"[PyramidRebuild] downsample_factors: {oz.metadata.downsample_factors}")
    print(f"[PyramidRebuild] downsample_method: {oz.metadata.downsample_method}")
    print(f"[PyramidRebuild] level 0 shape: {oz.shape(0)}")

    oz.rebuild_pyramid(n_threads=args.n_threads)

    print("[PyramidRebuild] Pyramid rebuild complete.")

if __name__ == "__main__":
    main()
