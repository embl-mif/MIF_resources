#!/usr/bin/env python3
"""Ultrack cell tracking wrapper, zarr v2 only compatible: consumes either an ilastik probability
store (foreground+contours) or Cellpose/SAM instance labels, tracks the
full timelapse in one job, and writes tracks.csv + a lineage graph. Writes
segments directly into an already-preallocated OME-Zarr level 0 (see
PREALLOCATE_ZARR / POPULATE_PYRAMID modules).
API reference: https://royerlab.github.io/ultrack/api.html
"""
import sys
sys.modules["gurobipy"] = None

import argparse
from pathlib import Path

import numpy as np
import zarr

from ultrack import MainConfig, load_config, track, to_tracks_layer, tracks_to_zarr
from ultrack.imgproc import robust_invert


def _open_array(path: str, level: str = "0"):
    root = zarr.open(path, mode="r")
    return root[level] if isinstance(root, zarr.hierarchy.Group) and level in root else root


def _squeeze_singleton_z(arr: np.ndarray) -> np.ndarray:
    """If arr is (T, Z, Y, X) with Z==1, squeeze to (T, Y, X). Leaves
    3D (real z-depth) or already-2D arrays unchanged."""
    if arr.ndim == 4 and arr.shape[1] == 1:
        return arr.squeeze(axis=1)
    return arr


def _load_time_sliced(path: str, channel: int, t_min: int, t_max: int):
    """Open a zarr array, select a channel, and slice to [t_min, t_max), only the requested timepoints are read
    from disk."""
    arr = _open_array(path)
    if channel is not None:
        arr = arr[:, channel]
    shape_t = arr.shape[0]
    t_start = t_min if t_min is not None else 0
    t_end = t_max if t_max is not None else shape_t
    if not (0 <= t_start < t_end <= shape_t):
        raise ValueError(f"Bad t range [{t_start}, {t_end}) for extent {shape_t}.")
    return np.asarray(arr[t_start:t_end])


def _write_segments_to_preallocated(segments: np.ndarray, output_path: str):
    """Write ultrack's raw label array into level 0 of an already-preallocated
    OME-Zarr store (see PREALLOCATE_ZARR). Reshapes to 5D tczyx to match the
    store's convention, then writes the whole array in one shot."""
    ndim = segments.ndim
    if ndim == 3:
        segments = segments.reshape(segments.shape[0], 1, 1, *segments.shape[1:])
    elif ndim == 4:
        segments = segments.reshape(segments.shape[0], 1, *segments.shape[1:])
    else:
        raise ValueError(f"Unexpected segments array ndim={ndim}, shape={segments.shape}")

    out_root = zarr.open(output_path, mode="r+")
    out_store = out_root["0"] if isinstance(out_root, zarr.hierarchy.Group) and "0" in out_root else out_root

    if out_store.shape != segments.shape:
        raise ValueError(
            f"Preallocated store has shape {out_store.shape}, but ultrack produced "
            f"segments of shape {segments.shape}. Was PREALLOCATE_ZARR run with the "
            f"wrong num_channels or dimensions?"
        )

    out_store[:] = segments
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Run ultrack cell tracking on a timelapse.")

    parser.add_argument("--labels-zarr", type=str, default=None)
    parser.add_argument("--foreground-zarr", type=str, default=None)
    parser.add_argument("--contours-zarr", type=str, default=None)
    parser.add_argument("--foreground-channel", type=int, default=None)
    parser.add_argument("--contours-channel", type=int, default=None)
    parser.add_argument("--derive-contours-from-raw", type=str, default=None,
                         help="Path to the raw intensity OME-Zarr. If given (and "
                              "--contours-zarr is not), contours are computed via "
                              "ultrack.imgproc.robust_invert on this raw data.")
    parser.add_argument("--raw-channel", type=int, default=None,
                         help="Channel index to select from --derive-contours-from-raw, if multi-channel.")

    parser.add_argument("--config-toml", type=str, default=None)
    parser.add_argument("--working-dir", type=str, required=True)
    parser.add_argument("--n-threads", type=int, default=1, help="Thread count for segmentation/tracking.")
    parser.add_argument("--t-min", type=int, default=None,
                         help="First timepoint index to process (inclusive). Defaults to 0.")
    parser.add_argument("--t-max", type=int, default=None,
                         help="Last timepoint index to process (exclusive). Defaults to full extent.")
    parser.add_argument("--scale-z", type=float, default=1.0)
    parser.add_argument("--scale-y", type=float, default=1.0)
    parser.add_argument("--scale-x", type=float, default=1.0)
    parser.add_argument("--min-area", type=int, default=None,
                         help="Minimum segment/cell area (pixels). Discards segments smaller than this.")
    parser.add_argument("--max-area", type=int, default=None,
                         help="Maximum segment/cell area (pixels). Discards segments larger than this.")
    parser.add_argument("--max-distance", type=float, default=None,
                         help="Maximum linking distance between consecutive frames (same units as --scale-*).")
    parser.add_argument("--solution-gap", type=float, default=None,
                         help="ILP solver optimality gap (larger = faster/less memory, less optimal).")
    parser.add_argument("--time-limit", type=int, default=None,
                         help="ILP solver time limit in seconds.")
    parser.add_argument("--distance-weight", type=float, default=None,
                         help="Weight given to link distance in the ILP objective. "
                              "0 (default) ignores distance beyond the --max-distance cutoff; "
                              "higher values bias the solver toward shorter links when choosing "
                              "among ambiguous candidates.")

    parser.add_argument("--output-tracks-csv", type=str, required=True)
    parser.add_argument("--output-graph", type=str, required=True,
                         help=".npy path for the lineage/division graph -- maps each track id to its parent(s).")
    parser.add_argument("--output-segments-zarr", type=str, required=True,
                         help="Path to the segments OME-Zarr store. Must already exist -- see PREALLOCATE_ZARR.")
    parser.add_argument("--overwrite", choices=["all", "links", "solutions", "none"], default="none")

    args = parser.parse_args()

    if args.labels_zarr and (args.foreground_zarr or args.contours_zarr or args.derive_contours_from_raw):
        parser.error("--labels-zarr is mutually exclusive with --foreground-zarr/--contours-zarr/--derive-contours-from-raw.")
    if not args.labels_zarr:
        if not args.foreground_zarr:
            parser.error("Provide either --labels-zarr, or --foreground-zarr with contours "
                          "(--contours-zarr, or --derive-contours-from-raw).")
        if args.contours_zarr and args.derive_contours_from_raw:
            parser.error("--contours-zarr and --derive-contours-from-raw are mutually exclusive.")
        if not args.contours_zarr and not args.derive_contours_from_raw:
            parser.error("Provide either --contours-zarr or --derive-contours-from-raw alongside --foreground-zarr.")

    config = load_config(args.config_toml) if args.config_toml else MainConfig()
    config.data_config.working_dir = Path(args.working_dir)
    config.segmentation_config.n_workers = args.n_threads
    config.tracking_config.n_threads = args.n_threads
    config.linking_config.n_workers = args.n_threads

    if args.min_area is not None:
        config.segmentation_config.min_area = args.min_area
    if args.max_area is not None:
        config.segmentation_config.max_area = args.max_area
    if args.max_distance is not None:
        config.linking_config.max_distance = args.max_distance
    if args.distance_weight is not None:
        config.linking_config.distance_weight = args.distance_weight
    if args.solution_gap is not None:
        config.tracking_config.solution_gap = args.solution_gap
    if args.time_limit is not None:
        config.tracking_config.time_limit = args.time_limit

    scale = [args.scale_z, args.scale_y, args.scale_x]

    if args.labels_zarr:
        print(f"[Ultrack] Loading labels from {args.labels_zarr} (t=[{args.t_min},{args.t_max}))")
        labels = _load_time_sliced(args.labels_zarr, args.foreground_channel, args.t_min, args.t_max)
        labels = _squeeze_singleton_z(labels)
        print(f"[Ultrack] labels shape: {labels.shape}")
        track(config, labels=labels, scale=scale, overwrite=args.overwrite)
    else:
        print(f"[Ultrack] Loading foreground from {args.foreground_zarr} (t=[{args.t_min},{args.t_max}))")
        foreground = _load_time_sliced(args.foreground_zarr, args.foreground_channel, args.t_min, args.t_max)
        foreground = _squeeze_singleton_z(foreground)

        if args.contours_zarr:
            print(f"[Ultrack] Loading contours from {args.contours_zarr} (t=[{args.t_min},{args.t_max}))")
            contours = _load_time_sliced(args.contours_zarr, args.contours_channel, args.t_min, args.t_max)
            contours = _squeeze_singleton_z(contours)
        else:
            print(f"[Ultrack] Deriving contours via robust_invert on {args.derive_contours_from_raw} (t=[{args.t_min},{args.t_max}))")
            raw = _load_time_sliced(args.derive_contours_from_raw, args.raw_channel, args.t_min, args.t_max)
            raw = _squeeze_singleton_z(raw)
            print(f"[Ultrack] raw shape going into robust_invert: {raw.shape}, dtype: {raw.dtype}")
            has_z = raw.ndim == 4
            voxel_size_per_frame = scale if has_z else scale[1:]

            contours_frames = []
            for t in range(raw.shape[0]):
                frame_contours = robust_invert(raw[t], voxel_size=voxel_size_per_frame)
                contours_frames.append(frame_contours)
            contours = np.stack(contours_frames, axis=0)

        print(f"[Ultrack] foreground shape: {foreground.shape}, contours shape: {contours.shape}")
        track(config, foreground=foreground, contours=contours, scale=scale, overwrite=args.overwrite)

    print("[Ultrack] Tracking complete. Exporting tracks.csv")
    tracks_df, graph = to_tracks_layer(config)
    tracks_df.to_csv(args.output_tracks_csv, index=False)
    print(f"[Ultrack] Saving lineage graph to {args.output_graph}")
    np.save(args.output_graph, graph, allow_pickle=True)

    print("[Ultrack] Writing segments into preallocated store")
    segments = tracks_to_zarr(config, tracks_df, overwrite=True)
    _write_segments_to_preallocated(np.asarray(segments), args.output_segments_zarr)

    print(f"[Ultrack] Done. Tracks: {args.output_tracks_csv}, Segments: {args.output_segments_zarr}")


if __name__ == "__main__":
    main()
