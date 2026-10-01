#!/usr/bin/env python3

"""
Otsu-threshold one channel of a PyMIF OME-Zarr at a chosen resolution level
and label the thresholded mask with skimage.measure.label, then rebuild the
full multiscale pyramid (same number of levels, chunksize and effective
downsample factors as the input) by nearest-neighbor resizing that canonical
labeling to every other level's exact shape -- this preserves label IDs
across levels instead of re-segmenting (and therefore re-labeling
inconsistently) at each resolution independently. Writes the result as a new
label OME-Zarr (no channel axis, data_type="label").
"""

import argparse
import copy
from pathlib import Path

import dask.array as da
import numpy as np
import pymif.microscope_manager as mm
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.transform import resize


def parse_args():
    parser = argparse.ArgumentParser(
        description="Otsu-threshold + label one channel of a PyMIF OME-Zarr."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--resolution-level",
        type=int,
        default=0,
        help="Pyramid level to threshold+label at; every other level is derived from it.",
    )
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


def drop_axis(per_level_values, axis_index: int) -> list[tuple]:
    return [
        tuple(value for i, value in enumerate(level_values) if i != axis_index)
        for level_values in per_level_values
    ]


def resize_label_image(label_image: np.ndarray, target_shape: tuple[int, ...], t_index_out: int | None) -> np.ndarray:
    """Nearest-neighbor resize a label image to `target_shape`, per timepoint."""
    if label_image.shape == target_shape:
        return label_image

    if t_index_out is None:
        return resize(
            label_image, target_shape, order=0, preserve_range=True, anti_aliasing=False
        ).astype(label_image.dtype)

    resized = np.zeros(target_shape, dtype=label_image.dtype)
    frame_target_shape = tuple(size for axis, size in enumerate(target_shape) if axis != t_index_out)
    for t in range(label_image.shape[t_index_out]):
        indexer_src = [slice(None)] * label_image.ndim
        indexer_src[t_index_out] = t
        indexer_dst = [slice(None)] * len(target_shape)
        indexer_dst[t_index_out] = t

        frame = label_image[tuple(indexer_src)]
        resized[tuple(indexer_dst)] = resize(
            frame, frame_target_shape, order=0, preserve_range=True, anti_aliasing=False
        ).astype(label_image.dtype)

    return resized


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
    t_index_out = axes_no_c.index("t") if t_index is not None else None

    canonical_shape = tuple(size for axis, size in enumerate(level_data.shape) if axis != c_index)
    if t_index_out is not None:
        canonical_label_image = np.zeros(canonical_shape, dtype=label_dtype)
        for t, frame in zip(timepoints, label_frames):
            indexer = [slice(None)] * canonical_label_image.ndim
            indexer[t_index_out] = t
            canonical_label_image[tuple(indexer)] = frame
    else:
        canonical_label_image = label_frames[0]

    # Rebuild every pyramid level from the canonical labeling by
    # nearest-neighbor resizing to that level's exact (channel-dropped)
    # shape, so label IDs stay consistent across the whole pyramid instead
    # of being independently (and inconsistently) re-segmented per level.
    output_levels = []
    for level_index, level in enumerate(source.data):
        target_shape = tuple(size for axis, size in enumerate(level.shape) if axis != c_index)
        target_chunks = tuple(size for axis, size in enumerate(level.chunksize) if axis != c_index)

        if level_index == args.resolution_level:
            level_label_image = canonical_label_image
        else:
            level_label_image = resize_label_image(canonical_label_image, target_shape, t_index_out)

        output_levels.append(da.from_array(level_label_image, chunks=target_chunks))

    metadata["axes"] = axes_no_c
    metadata["size"] = drop_axis(metadata["size"], c_index)
    metadata["chunksize"] = drop_axis(metadata["chunksize"], c_index)
    # Scales/units are only tracked for spatial axes (never "c"), and every
    # pyramid level is kept, so they carry over unchanged.
    metadata.pop("channel_names", None)
    metadata.pop("channel_colors", None)
    metadata["dtype"] = str(label_dtype)
    metadata["data_type"] = "label"
    metadata["is_label"] = True

    output = mm.ArrayManager(output_levels, metadata)
    output.to_zarr(str(args.output), overwrite=True, drop_singleton=False)

    source.close()

    print(f"Otsu multiplier applied: {args.otsu_multiplier}")
    print(f"Labeled channel {args.channel} at resolution level {args.resolution_level}")
    print(f"Wrote {len(output_levels)}-level label OME-Zarr to {args.output}")


if __name__ == "__main__":
    main()
