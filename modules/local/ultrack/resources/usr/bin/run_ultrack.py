#!/usr/bin/env python3
"""Ultrack cell tracking wrapper: consumes either an ilastik probability
store (foreground+contours) or Cellpose/SAM instance labels, tracks the
full timelapse in one job, and writes tracks.csv + a segments OME-Zarr
skeleton (level 0 only, via pymif's ZarrManager).
API reference: https://royerlab.github.io/ultrack/api.html
"""

import argparse
from pathlib import Path

import numpy as np
import zarr
import pymif.microscope_manager as mm
from ultrack import MainConfig, load_config, track, to_tracks_layer, tracks_to_zarr


def _open_array(path: str, level: str = "0"):
    root = zarr.open(path, mode="r")
    return root[level] if isinstance(root, zarr.hierarchy.Group) and level in root else root


def _write_segments_skeleton(segments: np.ndarray, output_path: str, scale_zyx):
    """Write ultrack's raw label array as OME-Zarr level 0 only, via pymif's
    ZarrManager"""
    ndim = segments.ndim
    if ndim == 3:
        segments = segments.reshape(segments.shape[0], 1, 1, *segments.shape[1:])
    elif ndim == 4:
        segments = segments.reshape(segments.shape[0], 1, *segments.shape[1:])
    else:
        raise ValueError(f"Unexpected segments array ndim={ndim}, shape={segments.shape}")

    chunks = (1, 1, min(64, segments.shape[2]), min(256, segments.shape[3]), min(256, segments.shape[4]))

    z = mm.ZarrManager(
        output_path,
        mode="a",
        metadata={
            "size": [segments.shape],
            "chunksize": [chunks],
            "scales": [tuple([1.0, 1.0] + list(scale_zyx))],
            "units": ("second", "", "micrometer", "micrometer", "micrometer"),
            "axes": "tczyx",
            "dtype": str(segments.dtype),
            "data_type": "label",
        },
    )
    z.write_image_region(
        segments,
        t=slice(0, segments.shape[0]), c=slice(0, 1),
        z=slice(0, segments.shape[2]), y=slice(0, segments.shape[3]), x=slice(0, segments.shape[4]),
        level=0,
    )
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Run ultrack cell tracking on a timelapse.")

    parser.add_argument("--labels-zarr", type=str, default=None)
    parser.add_argument("--foreground-zarr", type=str, default=None)
    parser.add_argument("--contours-zarr", type=str, default=None)
    parser.add_argument("--foreground-channel", type=int, default=None)
    parser.add_argument("--contours-channel", type=int, default=None)

    parser.add_argument("--config-toml", type=str, default=None)
    parser.add_argument("--working-dir", type=str, required=True)

    parser.add_argument("--scale-z", type=float, default=1.0)
    parser.add_argument("--scale-y", type=float, default=1.0)
    parser.add_argument("--scale-x", type=float, default=1.0)

    parser.add_argument("--output-tracks-csv", type=str, required=True)
    parser.add_argument("--output-segments-zarr", type=str, required=True,
                         help="Level-0-only OME-Zarr; rebuild the pyramid downstream via REBUILD_PYRAMID.")
    parser.add_argument("--overwrite", choices=["all", "links", "solutions", "none"], default="none")

    args = parser.parse_args()

    if args.labels_zarr and (args.foreground_zarr or args.contours_zarr):
        parser.error("--labels-zarr is mutually exclusive with --foreground-zarr/--contours-zarr.")
    if bool(args.foreground_zarr) != bool(args.contours_zarr):
        parser.error("--foreground-zarr and --contours-zarr must be given together.")
    if not args.labels_zarr and not args.foreground_zarr:
        parser.error("Provide either --labels-zarr, or both --foreground-zarr and --contours-zarr.")

    config = load_config(args.config_toml) if args.config_toml else MainConfig()
    config.data_config.working_dir = Path(args.working_dir)

    scale = [args.scale_z, args.scale_y, args.scale_x]

    if args.labels_zarr:
        print(f"[Ultrack] Loading labels from {args.labels_zarr}")
        labels = _open_array(args.labels_zarr)
        if args.foreground_channel is not None:
            labels = labels[:, args.foreground_channel]
        print(f"[Ultrack] labels shape: {labels.shape}")
        track(config, labels=labels, scale=scale, overwrite=args.overwrite)
    else:
        print(f"[Ultrack] Loading foreground from {args.foreground_zarr}")
        foreground = _open_array(args.foreground_zarr)
        if args.foreground_channel is not None:
            foreground = foreground[:, args.foreground_channel]

        print(f"[Ultrack] Loading contours from {args.contours_zarr}")
        contours = _open_array(args.contours_zarr)
        if args.contours_channel is not None:
            contours = contours[:, args.contours_channel]

        print(f"[Ultrack] foreground shape: {foreground.shape}, contours shape: {contours.shape}")
        track(config, foreground=foreground, contours=contours, scale=scale, overwrite=args.overwrite)

    print("[Ultrack] Tracking complete. Exporting tracks.csv")
    tracks_df, graph = to_tracks_layer(config)
    tracks_df.to_csv(args.output_tracks_csv, index=False)

    print("[Ultrack] Writing segments skeleton (level 0 only, via pymif)")
    segments = tracks_to_zarr(config, tracks_df, overwrite=True)
    _write_segments_skeleton(np.asarray(segments), args.output_segments_zarr, scale)

    print(f"[Ultrack] Done. Tracks: {args.output_tracks_csv}, Segments: {args.output_segments_zarr}")


if __name__ == "__main__":
    main()
