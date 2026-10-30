#!/usr/bin/env python3

"""
Mask one channel of a ZarrSmith intensity OME-Zarr with an existing (already
filtered/relabeled) ZarrSmith label OME-Zarr -- zeroing every pixel outside a
labeled object -- and append the result as a new "clean" channel of a new
OME-Zarr holding every original channel plus this one. The intensity and
label stores must share the same multiscale pyramid (same number of levels
and matching shapes level-for-level, as produced by zarrsmith_otsu_label +
zarrsmith_filter); every level is masked and written out, so the output keeps
the same pyramid structure as the input intensity zarr.
"""

import argparse
import copy
from pathlib import Path

import dask.array as da
import numpy as np
import zarrsmith.microscope_manager as mm


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Zero one intensity channel outside the objects of a label "
            "OME-Zarr, and append it as a new channel of a new OME-Zarr."
        )
    )
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--channel", required=True, type=int)
    return parser.parse_args()


def extract_frame(level_data, indices: dict[int, int | None]) -> np.ndarray:
    indexer = [slice(None)] * level_data.ndim
    for axis, value in indices.items():
        if value is not None:
            indexer[axis] = value
    return np.asarray(level_data[tuple(indexer)])


def match_shape(array: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    """Crop/zero-pad `array` at the high end of each axis to `target_shape`.

    Independently-built multiscale pyramids (e.g. ZarrSmith's own converter vs.
    ngff-zarr's to_multiscales, used by Cellpose) can round a coarse level's
    dimensions a pixel or two differently even from the same base image.
    Mirrors _match_mask_shape in cellpose_sam_zarr.py.
    """
    if array.shape == target_shape:
        return array
    if array.ndim != len(target_shape):
        raise ValueError(f"Rank mismatch: array shape {array.shape} vs target {target_shape}")

    adjusted = array
    for axis, target in enumerate(target_shape):
        current = adjusted.shape[axis]
        if current > target:
            slicer = [slice(None)] * adjusted.ndim
            slicer[axis] = slice(0, target)
            adjusted = adjusted[tuple(slicer)]
        elif current < target:
            pad_width = [(0, 0)] * adjusted.ndim
            pad_width[axis] = (0, target - current)
            adjusted = np.pad(adjusted, pad_width=pad_width, mode="constant", constant_values=0)

    return adjusted


def mask_level(
    intensity_level,
    label_level,
    c_index: int,
    channel_idx: int,
    t_index_int: int | None,
    t_index_lbl: int | None,
) -> np.ndarray:
    """Zero `channel_idx` of `intensity_level` outside `label_level`'s objects.

    Returns a full-rank array matching `intensity_level`'s shape, with a
    size-1 channel axis at `c_index`.
    """
    masked_channel = np.zeros(
        tuple(
            1 if axis == c_index else intensity_level.shape[axis]
            for axis in range(intensity_level.ndim)
        ),
        dtype=intensity_level.dtype,
    )

    timepoints = range(intensity_level.shape[t_index_int]) if t_index_int is not None else [None]
    for t in timepoints:
        channel_image = extract_frame(
            intensity_level,
            {c_index: channel_idx, **({t_index_int: t} if t_index_int is not None else {})},
        )
        label_frame = extract_frame(
            label_level,
            {t_index_lbl: t} if t_index_lbl is not None else {},
        )

        if channel_image.shape != label_frame.shape:
            print(
                f"Adjusting label frame shape from {label_frame.shape} to "
                f"{channel_image.shape} (pyramid rounding mismatch between "
                "--intensity and --labels)"
            )
            label_frame = match_shape(label_frame, channel_image.shape)

        masked = np.where(label_frame != 0, channel_image, 0).astype(channel_image.dtype)

        indexer = [slice(None)] * masked_channel.ndim
        indexer[c_index] = 0
        if t_index_int is not None:
            indexer[t_index_int] = t
        masked_channel[tuple(indexer)] = masked

    return masked_channel


def main():
    args = parse_args()

    if args.output.exists():
        raise FileExistsError(f"Output already exists: {args.output}")

    intensity_source = mm.ZarrManager(str(args.intensity))
    metadata = copy.deepcopy(intensity_source.metadata)
    axes = metadata["axes"].lower()

    if "c" not in axes:
        raise ValueError(f"--intensity dataset does not contain a channel axis: axes={axes}")
    c_index = axes.index("c")
    n_channels = intensity_source.data[0].shape[c_index]
    if args.channel < 0 or args.channel >= n_channels:
        raise ValueError(f"--channel {args.channel} is out of range for {n_channels} channel(s)")

    labels_source = mm.ZarrManager(str(args.labels))
    labels_axes = labels_source.metadata["axes"].lower()

    if len(intensity_source.data) != len(labels_source.data):
        raise ValueError(
            f"--intensity has {len(intensity_source.data)} pyramid level(s) but "
            f"--labels has {len(labels_source.data)}. They must share the same "
            "multiscale pyramid (e.g. --labels produced by zarrsmith_otsu_label + "
            "zarrsmith_filter from this same --intensity zarr)."
        )

    t_index_int = axes.index("t") if "t" in axes else None
    t_index_lbl = labels_axes.index("t") if "t" in labels_axes else None
    if (t_index_int is not None) != (t_index_lbl is not None):
        raise ValueError(
            "--intensity and --labels disagree on whether a time axis is "
            f"present (intensity axes={axes}, labels axes={labels_axes})"
        )

    output_levels = []
    for intensity_level, label_level in zip(intensity_source.data, labels_source.data):
        masked_channel = mask_level(
            intensity_level, label_level, c_index, args.channel, t_index_int, t_index_lbl
        )
        target_chunks = tuple(
            1 if axis == c_index else size
            for axis, size in enumerate(intensity_level.chunksize)
        )
        combined = da.concatenate(
            [intensity_level, da.from_array(masked_channel, chunks=target_chunks)],
            axis=c_index,
        )
        output_levels.append(combined)

    clean_channel_name = "channel"
    channel_names = metadata.get("channel_names")
    if channel_names and args.channel < len(channel_names):
        clean_channel_name = str(channel_names[args.channel])
    new_channel_name = f"{clean_channel_name}_clean"

    if channel_names:
        metadata["channel_names"] = list(channel_names) + [new_channel_name]

    channel_colors = metadata.get("channel_colors")
    if channel_colors:
        source_color = (
            channel_colors[args.channel] if args.channel < len(channel_colors) else "#FFFFFF"
        )
        metadata["channel_colors"] = list(channel_colors) + [source_color]

    metadata["size"] = [tuple(int(size) for size in level.shape) for level in output_levels]
    metadata["chunksize"] = [tuple(int(size) for size in level.chunksize) for level in output_levels]
    metadata["dtype"] = str(output_levels[0].dtype)
    metadata["data_type"] = "intensity"

    output = mm.ArrayManager(output_levels, metadata)
    output.to_zarr(str(args.output), overwrite=True, drop_singleton=False)

    intensity_source.close()
    labels_source.close()

    print(f"New clean channel appended at index {n_channels}: '{new_channel_name}'")
    print(f"Wrote {len(output_levels)}-level OME-Zarr to {args.output}")


if __name__ == "__main__":
    main()
