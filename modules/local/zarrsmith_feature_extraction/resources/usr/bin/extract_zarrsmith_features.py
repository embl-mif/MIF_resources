#!/usr/bin/env python3

"""
Load an intensity OME-Zarr and a label OME-Zarr with ZarrSmith, measure every
label with skimage.measure.regionprops (a base property set, plus any
user-requested structural and/or per-channel intensity properties), and
write the resulting feature table as a single CSV.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import zarrsmith.microscope_manager as mm
from skimage.measure import regionprops_table

# Shape/geometry properties always measured (when supported for the label
# image's dimensionality). "perimeter"/"orientation" are 2D-only in skimage
# and are used to derive "form_factor" (4*pi*area/perimeter^2).
BASE_STRUCTURAL_PROPERTIES_2D = (
    "label", "area", "bbox", "centroid", "perimeter", "orientation",
    "eccentricity", "solidity", "extent",
    "major_axis_length", "minor_axis_length",
)
BASE_STRUCTURAL_PROPERTIES_ND = (
    "label", "area", "bbox", "centroid", "solidity", "extent",
)

BASE_INTENSITY_PROPERTIES = ("label", "mean_intensity", "min_intensity", "max_intensity")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Measure skimage regionprops for every label of a label OME-Zarr, "
            "against every channel of an intensity OME-Zarr, and write the "
            "features to a CSV."
        )
    )
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--labels-resolution-level", type=int, default=0)
    parser.add_argument(
        "--channels",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated intensity channel indices "
            "to measure. Defaults to every channel."
        ),
    )
    parser.add_argument(
        "--extra-properties",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated extra skimage.measure."
            "regionprops properties to add to the base structural set "
            "(measured once per label, not per channel)."
        ),
    )
    parser.add_argument(
        "--extra-intensity-properties",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated extra skimage.measure."
            "regionprops properties that require an intensity image, "
            "measured once per requested channel (e.g. 'weighted_centroid')."
        ),
    )
    parser.add_argument("--dataset-id", type=str, default=None)
    return parser.parse_args()


def parse_token_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [token for token in re.split(r"[;,\s]+", text.strip()) if token]


def parse_int_list(text: str | None) -> list[int]:
    tokens = parse_token_list(text)
    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid channel list '{text}'. Expected integer channel "
            "indices separated by comma, semicolon, or whitespace."
        ) from exc


def sanitize_suffix(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_") or "channel"


def measure_structural_regions(label_image: np.ndarray, extra_properties: list[str]) -> pd.DataFrame:
    is_2d = label_image.ndim == 2
    base_properties = BASE_STRUCTURAL_PROPERTIES_2D if is_2d else BASE_STRUCTURAL_PROPERTIES_ND

    try:
        table = regionprops_table(label_image, properties=base_properties)
    except NotImplementedError:
        base_properties = BASE_STRUCTURAL_PROPERTIES_ND
        table = regionprops_table(label_image, properties=base_properties)

    df = pd.DataFrame(table)

    if "perimeter" in df.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["form_factor"] = np.where(
                df["perimeter"] > 0,
                4 * np.pi * df["area"] / np.square(df["perimeter"]),
                0.0,
            )

    if extra_properties:
        # Requested explicitly by the user: let a NotImplementedError (e.g. a
        # 2D-only property against a 3D label image) surface as a clear error
        # rather than being silently dropped like the opportunistic base set.
        extra_table = regionprops_table(label_image, properties=("label", *extra_properties))
        df = df.merge(pd.DataFrame(extra_table), on="label")

    return df


def measure_intensity(
    label_image: np.ndarray,
    intensity_image: np.ndarray,
    suffix: str,
    extra_intensity_properties: list[str],
) -> pd.DataFrame:
    properties = BASE_INTENSITY_PROPERTIES + tuple(extra_intensity_properties)
    table = regionprops_table(label_image, intensity_image=intensity_image, properties=properties)
    df = pd.DataFrame(table)
    return df.rename(columns={col: f"{col}_{suffix}" for col in df.columns if col != "label"})


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


def main():
    args = parse_args()

    channel_indices = parse_int_list(args.channels)
    extra_properties = parse_token_list(args.extra_properties)
    extra_intensity_properties = parse_token_list(args.extra_intensity_properties)

    intensity_source = mm.ZarrManager(str(args.intensity))
    intensity_axes = intensity_source.metadata["axes"].lower()
    if args.intensity_resolution_level < 0 or args.intensity_resolution_level >= len(intensity_source.data):
        raise ValueError(
            f"--intensity-resolution-level {args.intensity_resolution_level} is out of "
            f"range for input with {len(intensity_source.data)} level(s)"
        )
    intensity_data = intensity_source.data[args.intensity_resolution_level]

    if "c" not in intensity_axes:
        raise ValueError(f"Intensity dataset does not contain a channel axis: axes={intensity_axes}")
    c_index = intensity_axes.index("c")
    n_channels = intensity_data.shape[c_index]

    if not channel_indices:
        channel_indices = list(range(n_channels))
    for channel_idx in channel_indices:
        if channel_idx < 0 or channel_idx >= n_channels:
            raise ValueError(f"Channel index {channel_idx} is out of range for {n_channels} channel(s)")

    channel_names = intensity_source.metadata.get("channel_names")
    channel_suffixes = {}
    for channel_idx in channel_indices:
        if channel_names and channel_idx < len(channel_names):
            channel_suffixes[channel_idx] = f"ch{channel_idx}_{sanitize_suffix(str(channel_names[channel_idx]))}"
        else:
            channel_suffixes[channel_idx] = f"ch{channel_idx}"

    labels_source = mm.ZarrManager(str(args.labels))
    labels_axes = labels_source.metadata["axes"].lower()
    if args.labels_resolution_level < 0 or args.labels_resolution_level >= len(labels_source.data):
        raise ValueError(
            f"--labels-resolution-level {args.labels_resolution_level} is out of "
            f"range for input with {len(labels_source.data)} level(s)"
        )
    labels_data = labels_source.data[args.labels_resolution_level]

    t_index_int = intensity_axes.index("t") if "t" in intensity_axes else None
    t_index_lbl = labels_axes.index("t") if "t" in labels_axes else None
    # A singleton T axis on the intensity image is dropped when labels have none.
    squeeze_int_t = (
        t_index_int is not None
        and t_index_lbl is None
        and intensity_data.shape[t_index_int] == 1
    )
    if squeeze_int_t:
        print(f"Dropping singleton T axis from intensity image (axes={intensity_axes})")
    elif (t_index_int is not None) != (t_index_lbl is not None):
        raise ValueError(
            "Intensity and labels disagree on whether a time axis is present "
            f"(intensity axes={intensity_axes}, labels axes={labels_axes})"
        )
    if (
        not squeeze_int_t
        and t_index_int is not None
        and intensity_data.shape[t_index_int] != labels_data.shape[t_index_lbl]
    ):
        raise ValueError(
            "Intensity and labels have a different number of timepoints "
            f"({intensity_data.shape[t_index_int]} vs {labels_data.shape[t_index_lbl]})"
        )

    if squeeze_int_t:
        timepoints = [None]
        int_t_fixed = {t_index_int: 0}
    else:
        timepoints = range(intensity_data.shape[t_index_int]) if t_index_int is not None else [None]
        int_t_fixed = {}

    frame_tables = []

    for t in timepoints:
        label_image = extract_frame(labels_data, {t_index_lbl: t} if t_index_lbl is not None else {})
        label_image = label_image.astype(np.int64, copy=False)

        first_channel_image = extract_frame(
            intensity_data,
            {c_index: channel_indices[0], **int_t_fixed, **({t_index_int: t} if t is not None and t_index_int is not None else {})},
        )
        if label_image.shape != first_channel_image.shape:
            print(
                f"Adjusting label image shape from {label_image.shape} to "
                f"{first_channel_image.shape} (pyramid rounding mismatch between "
                f"--intensity-resolution-level={args.intensity_resolution_level} and "
                f"--labels-resolution-level={args.labels_resolution_level})"
            )
            label_image = match_shape(label_image, first_channel_image.shape)

        df = measure_structural_regions(label_image, extra_properties)

        for channel_idx in channel_indices:
            channel_image = (
                first_channel_image
                if channel_idx == channel_indices[0]
                else extract_frame(
                    intensity_data,
                    {c_index: channel_idx, **int_t_fixed, **({t_index_int: t} if t is not None and t_index_int is not None else {})},
                )
            )
            intensity_df = measure_intensity(
                label_image, channel_image, channel_suffixes[channel_idx], extra_intensity_properties
            )
            df = df.merge(intensity_df, on="label")

        if t is not None:
            df.insert(0, "timepoint", t)
        frame_tables.append(df)

    features = pd.concat(frame_tables, ignore_index=True) if frame_tables else pd.DataFrame()
    if args.dataset_id:
        features.insert(0, "dataset_id", args.dataset_id)

    features.to_csv(args.output, index=False)

    intensity_source.close()
    labels_source.close()

    print(f"Measured channel(s): {[(idx, channel_suffixes[idx]) for idx in channel_indices]}")
    print(f"Wrote {len(features)} row(s) across {len(timepoints)} timepoint(s) to {args.output}")


if __name__ == "__main__":
    main()
