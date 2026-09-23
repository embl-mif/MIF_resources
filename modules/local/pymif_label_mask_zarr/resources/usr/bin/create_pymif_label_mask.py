#!/usr/bin/env python3

"""
Threshold one channel of a PyMIF OME-Zarr with Otsu, label the thresholded
mask with skimage, filter the resulting objects by area, fluorescence
intensity in any channel, and/or any other skimage regionprops property, use
the filtered label mask to zero out the source channel outside the kept
objects, and write a new single-resolution OME-Zarr containing every original
channel plus this new masked channel appended at the end.
"""

import argparse
import copy
import re
import warnings
from pathlib import Path
from typing import Any

import dask.array as da
import numpy as np
import pandas as pd
import pymif.microscope_manager as mm
from skimage.measure import label, regionprops_table
from skimage.filters import threshold_otsu

# Shape/geometry regionprops properties always measured (when supported for
# the image's dimensionality). "perimeter" is 2D-only in skimage and is used
# to derive "form_factor" (4*pi*area/perimeter^2, a roundness measure in
# [0, 1]). Per-channel intensity properties (mean/min/max_intensity) are
# measured separately, once per requested channel -- see measure_intensity().
STRUCTURAL_PROPERTIES_2D = (
    "label", "area", "bbox", "centroid", "perimeter",
    "eccentricity", "solidity", "extent",
    "major_axis_length", "minor_axis_length",
)
STRUCTURAL_PROPERTIES_ND = (
    "label", "area", "bbox", "centroid", "solidity", "extent",
)

INTENSITY_PROPERTIES = ("label", "mean_intensity", "min_intensity", "max_intensity")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Otsu-threshold + label + filter one channel, mask it, and "
            "append it as a new channel of a new OME-Zarr."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resolution-level", type=int, default=0)
    parser.add_argument("--mask-channel", required=True, type=int)
    parser.add_argument(
        "--otsu-multiplier",
        type=float,
        default=1.0,
        help="Multiplier applied to the Otsu threshold before binarizing.",
    )
    parser.add_argument("--min-area", type=float, default=None)
    parser.add_argument("--max-area", type=float, default=None)
    parser.add_argument(
        "--intensity-channels",
        type=str,
        default=None,
        help=(
            "Comma/semicolon/whitespace-separated channel indices to also "
            "measure per-object fluorescence intensity on (in addition to "
            "--mask-channel). Each channel c yields "
            "mean_intensity_ch<c>/min_intensity_ch<c>/max_intensity_ch<c> "
            "columns, filterable via --filters."
        ),
    )
    parser.add_argument(
        "--filters",
        type=str,
        default=None,
        help=(
            "Semicolon-separated 'property:min:max' bounds on any measured "
            "regionprops property (including 'form_factor' and the "
            "mean/min/max_intensity[_ch<c>] columns). Either bound may be "
            "left empty, e.g. 'solidity:0.8:;mean_intensity_ch1:100:'."
        ),
    )
    return parser.parse_args()


def parse_filters(text: str | None) -> dict[str, tuple[float | None, float | None]]:
    if not text:
        return {}

    filters: dict[str, tuple[float | None, float | None]] = {}
    for token in text.split(";"):
        token = token.strip()
        if not token:
            continue

        parts = token.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid filter '{token}'. Expected 'property:min:max' "
                "(either bound may be empty)."
            )

        prop, lo_text, hi_text = (part.strip() for part in parts)
        if not prop:
            raise ValueError(f"Invalid filter '{token}': missing property name.")

        lo = float(lo_text) if lo_text else None
        hi = float(hi_text) if hi_text else None
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(
                f"Invalid filter '{token}': max must be greater than min."
            )

        filters[prop] = (lo, hi)

    return filters


def parse_int_list(text: str | None) -> list[int]:
    if not text:
        return []

    tokens = [token for token in re.split(r"[;,\s]+", text.strip()) if token]
    try:
        return [int(token) for token in tokens]
    except ValueError as exc:
        raise ValueError(
            f"Invalid channel list '{text}'. Expected integer channel "
            "indices separated by comma, semicolon, or whitespace."
        ) from exc


def measure_structural_regions(label_image: np.ndarray) -> pd.DataFrame:
    is_2d = label_image.ndim == 2
    properties = STRUCTURAL_PROPERTIES_2D if is_2d else STRUCTURAL_PROPERTIES_ND

    try:
        table = regionprops_table(label_image, properties=properties)
    except NotImplementedError:
        warnings.warn(
            f"Some regionprops properties are unsupported for {label_image.ndim}D "
            "images; falling back to a reduced property set.",
            stacklevel=2,
        )
        properties = STRUCTURAL_PROPERTIES_ND
        table = regionprops_table(label_image, properties=properties)

    df = pd.DataFrame(table)

    if "perimeter" in df.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            df["form_factor"] = np.where(
                df["perimeter"] > 0,
                4 * np.pi * df["area"] / np.square(df["perimeter"]),
                0.0,
            )

    return df


def measure_intensity(label_image: np.ndarray, intensity_image: np.ndarray, channel_idx: int) -> pd.DataFrame:
    table = regionprops_table(
        label_image,
        intensity_image=intensity_image,
        properties=INTENSITY_PROPERTIES,
    )
    df = pd.DataFrame(table)
    return df.rename(columns={
        "mean_intensity": f"mean_intensity_ch{channel_idx}",
        "min_intensity": f"min_intensity_ch{channel_idx}",
        "max_intensity": f"max_intensity_ch{channel_idx}",
    })


def filter_regions(
    df: pd.DataFrame,
    min_area: float | None,
    max_area: float | None,
    filters: dict[str, tuple[float | None, float | None]],
) -> pd.DataFrame:
    keep = pd.Series(True, index=df.index)

    if min_area is not None:
        keep &= df["area"] >= min_area
    if max_area is not None:
        keep &= df["area"] <= max_area

    for prop, (lo, hi) in filters.items():
        if prop not in df.columns:
            raise ValueError(
                f"Unknown or unmeasured region property '{prop}'. Available "
                f"properties: {sorted(df.columns)}"
            )
        if lo is not None:
            keep &= df[prop] >= lo
        if hi is not None:
            keep &= df[prop] <= hi

    return df[keep]


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

    filters = parse_filters(args.filters)
    requested_intensity_channels = parse_int_list(args.intensity_channels)

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

    if args.mask_channel < 0 or args.mask_channel >= n_channels:
        raise ValueError(
            f"--mask-channel {args.mask_channel} is out of range for "
            f"{n_channels} channel(s)"
        )

    for channel_idx in requested_intensity_channels:
        if channel_idx < 0 or channel_idx >= n_channels:
            raise ValueError(
                f"--intensity-channels index {channel_idx} is out of range "
                f"for {n_channels} channel(s)"
            )

    # The mask channel's own intensity is always measured (as both the plain
    # mean/min/max_intensity columns and the mean/min/max_intensity_ch<c>
    # columns), in addition to any other channels requested for filtering.
    intensity_channels = sorted({args.mask_channel, *requested_intensity_channels})

    t_index = axes.index("t") if "t" in axes else None
    timepoints = range(level_data.shape[t_index]) if t_index is not None else [None]

    masked_frames = []
    kept_tables = []

    for t in timepoints:
        channel_image = extract_channel_image(level_data, c_index, t_index, args.mask_channel, t)

        threshold = threshold_otsu(channel_image) * args.otsu_multiplier
        binary = channel_image > threshold
        label_image = label(binary)

        df = measure_structural_regions(label_image)
        for channel_idx in intensity_channels:
            intensity_image = (
                channel_image
                if channel_idx == args.mask_channel
                else extract_channel_image(level_data, c_index, t_index, channel_idx, t)
            )
            df = df.merge(measure_intensity(label_image, intensity_image, channel_idx), on="label")

        # Plain aliases for the mask channel's own intensity, for convenience.
        for stat in ("mean_intensity", "min_intensity", "max_intensity"):
            df[stat] = df[f"{stat}_ch{args.mask_channel}"]

        kept = filter_regions(df, args.min_area, args.max_area, filters)

        keep_mask = np.isin(label_image, kept["label"].to_numpy())
        masked_image = np.where(keep_mask, channel_image, 0).astype(channel_image.dtype)

        if t is not None:
            kept = kept.copy()
            kept.insert(0, "timepoint", t)
        kept_tables.append(kept)
        masked_frames.append(masked_image)

    # Assemble the new channel back into a full-rank array (matching
    # level_data's shape with a size-1 channel axis) regardless of where the
    # "t" and "c" axes sit relative to each other and to the spatial axes.
    masked_channel = np.zeros(
        tuple(
            1 if axis == c_index else level_data.shape[axis]
            for axis in range(level_data.ndim)
        ),
        dtype=level_data.dtype,
    )
    for t, masked_image in zip(timepoints, masked_frames):
        indexer = [slice(None)] * level_data.ndim
        indexer[c_index] = 0
        if t_index is not None:
            indexer[t_index] = t
        masked_channel[tuple(indexer)] = masked_image

    combined = da.concatenate(
        [level_data, da.from_array(masked_channel, chunks=level_data.chunksize)],
        axis=c_index,
    )

    mask_channel_name = "channel"
    channel_names = metadata.get("channel_names")
    if channel_names and args.mask_channel < len(channel_names):
        mask_channel_name = str(channel_names[args.mask_channel])
    new_channel_name = f"{mask_channel_name}_masked"

    if channel_names:
        metadata["channel_names"] = list(channel_names) + [new_channel_name]

    channel_colors = metadata.get("channel_colors")
    if channel_colors:
        source_color = (
            channel_colors[args.mask_channel]
            if args.mask_channel < len(channel_colors)
            else "#FFFFFF"
        )
        metadata["channel_colors"] = list(channel_colors) + [source_color]

    metadata["size"] = [tuple(int(size) for size in combined.shape)]
    metadata["chunksize"] = [tuple(int(size) for size in combined.chunksize)]
    if "scales" in metadata:
        metadata["scales"] = [metadata["scales"][args.resolution_level]]
    metadata["dtype"] = str(combined.dtype)
    metadata["data_type"] = "intensity"

    output = mm.ArrayManager([combined], metadata)
    output.to_zarr(str(args.output), overwrite=True)

    source.close()

    kept_all = pd.concat(kept_tables, ignore_index=True) if kept_tables else pd.DataFrame()
    kept_all.insert(0, "mask_channel", mask_channel_name)
    measurements_path = args.output.with_name(f"{args.output.stem}_regionprops.csv")
    kept_all.to_csv(measurements_path, index=False)

    print(f"Otsu multiplier applied: {args.otsu_multiplier}")
    print(f"Measured intensity for channel(s): {intensity_channels}")
    print(f"New masked channel appended at index {n_channels}: '{new_channel_name}'")
    print(f"Kept {len(kept_all)} object(s) across {len(timepoints)} timepoint(s)")
    print(f"Wrote regionprops measurements to {measurements_path}")


if __name__ == "__main__":
    main()
