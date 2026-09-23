#!/usr/bin/env python3

"""
Otsu-threshold one channel of a PyMIF OME-Zarr and label the thresholded
mask with skimage.measure.label, and write the result as a new
single-resolution label OME-Zarr (no channel axis, data_type="label").
"""

import argparse
import copy
from pathlib import Path

import numpy as np
import pymif.microscope_manager as mm
from skimage.filters import threshold_otsu
from skimage.measure import label


def parse_args():
    parser = argparse.ArgumentParser(
        description="Otsu-threshold + label one channel of a PyMIF OME-Zarr."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resolution-level", type=int, default=0)
    parser.add_argument("--channel", required=True, type=int)
    parser.add_argument(
        "--otsu-multiplier",
        type=float,
        default=1.0,
        help="Multiplier applied to the Otsu threshold before binarizing.",
    )
    return parser.parse_args()


def extract_channel_image(level_data, c_index: int, t_index: int | None, channel_idx: int, t: int | None) -> np.ndarray:
    indexer = [slice(None)] * level_data.ndim
    indexer[c_index] = channel_idx
    if t_index is not None:
        indexer[t_index] = t
    return np.asarray(level_data[tuple(indexer)])


def main():
    args = parse_args()

    if args.output.exists():
        raise FileExistsError(f"Output already exists: {args.output}")

    source = mm.ZarrManager(str(args.input))
    metadata = copy.deepcopy(source.metadata)
    axes = metadata["axes"].lower()

    if args.resolution_level < 0 or args.resolution_level >= len(source.data):
        raise ValueError(
            f"--resolution-level {args.resolution_level} is out of range for "
            f"input with {len(source.data)} level(s)"
        )
    level_data = source.data[args.resolution_level]

    if "c" not in axes:
        raise ValueError(f"Input dataset does not contain a channel axis: axes={axes}")
    c_index = axes.index("c")
    n_channels = level_data.shape[c_index]

    if args.channel < 0 or args.channel >= n_channels:
        raise ValueError(f"--channel {args.channel} is out of range for {n_channels} channel(s)")

    t_index = axes.index("t") if "t" in axes else None
    timepoints = range(level_data.shape[t_index]) if t_index is not None else [None]

    label_frames = []
    for t in timepoints:
        channel_image = extract_channel_image(level_data, c_index, t_index, args.channel, t)
        threshold = threshold_otsu(channel_image) * args.otsu_multiplier
        binary = channel_image > threshold
        label_frames.append(label(binary))

    label_dtype = label_frames[0].dtype

    # Output axes drop the channel axis entirely -- label stores have no
    # channel dimension (matching e.g. the Cellpose-SAM label output in this
    # repo). Spatial/time axes and their relative order are unchanged.
    axes_no_c = axes.replace("c", "")
    if t_index is not None:
        t_index_out = axes_no_c.index("t")
        out_shape = tuple(
            level_data.shape[axis] if axis != c_index else None
            for axis in range(level_data.ndim)
        )
        out_shape = tuple(size for size in out_shape if size is not None)
        label_image = np.zeros(out_shape, dtype=label_dtype)
        for t, frame in zip(timepoints, label_frames):
            indexer = [slice(None)] * label_image.ndim
            indexer[t_index_out] = t
            label_image[tuple(indexer)] = frame
    else:
        label_image = label_frames[0]

    metadata["axes"] = axes_no_c
    metadata["size"] = [tuple(int(size) for size in label_image.shape)]
    metadata["chunksize"] = [tuple(int(size) for size in label_image.shape)]
    if "scales" in metadata:
        # Scales/units are only tracked for spatial axes (never "c"), so
        # dropping the channel axis leaves them unchanged; just select the
        # single resolution level being written out.
        metadata["scales"] = [metadata["scales"][args.resolution_level]]
    metadata.pop("channel_names", None)
    metadata.pop("channel_colors", None)
    metadata["dtype"] = str(label_dtype)
    metadata["data_type"] = "label"
    metadata["is_label"] = True

    output = mm.ArrayManager([label_image], metadata)
    output.to_zarr(str(args.output), overwrite=True)

    source.close()

    print(f"Otsu multiplier applied: {args.otsu_multiplier}")
    print(f"Labeled channel {args.channel} at resolution level {args.resolution_level}")
    print(f"Wrote label OME-Zarr to {args.output}")


if __name__ == "__main__":
    main()
