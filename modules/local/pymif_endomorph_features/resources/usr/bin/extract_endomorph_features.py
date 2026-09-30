#!/usr/bin/env python3

"""
Measure per-cell morphology and per-channel band intensities against the
regions OME-Zarr written by create_endomorph_regions.py, and summarize the
image, including its empty area.

With N ring depths (--ring-depths, only the count is used here), the regions
zarr holds bands 1..N (R1..RN, nested: RK = pixels with value 1..K), N+1 =
empty and N+2 = not imaged. A cell is measured if create_endomorph_regions.py
painted any band pixel for it; a cell with no painted band pixel was excluded
there for touching the image border.

Writes one CSV of per-cell measurements (morphology + the mean of every
--measure-channels channel in each of R1..RN and in the whole cell) and one
CSV of per-image summary metrics (cell counts and empty-area statistics).

Stripped down from endo_morph's features.py/run_endo_morph.py
(Z:\\people\\gritti\\projects\\Candelas_Adrian\\endo_morph).
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pymif.microscope_manager as mm
from skimage.measure import regionprops


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Measure per-cell morphology and band intensities against the "
            "R1..RN bands of --regions, and summarize the image (cell counts "
            "and empty area); write both as CSVs."
        )
    )
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--ring-depths", type=str, default="6 9 24",
                        help="The ring depths --regions was made with; only their number N "
                             "is used (bands 1..N, N+1 = empty, N+2 = not imaged).")
    parser.add_argument("--measure-channels", type=str, default=None,
                        help="Intensity channels (index or exact name, separated by "
                             "comma/semicolon/whitespace) to average in every region. "
                             "Default: all channels.")
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--labels-resolution-level", type=int, default=0)
    parser.add_argument("--timepoint", type=int, default=0)
    parser.add_argument("--dataset-id", type=str, default=None)
    return parser.parse_args()


def parse_token_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [token for token in re.split(r"[;,\s]+", text.strip()) if token]


# --- loading -----------------------------------------------------------------


def extract_plane(level_data, axes: str, t: int, keep_channel: bool) -> np.ndarray:
    """Select timepoint `t`, max-project Z if present, and return (C, Y, X) or (Y, X)."""
    axes = axes.lower()
    indexer = [t if ax == "t" else slice(None) for ax in axes]
    arr = level_data[tuple(indexer)]
    axes = axes.replace("t", "")

    if "z" in axes:
        z_ax = axes.index("z")
        arr = arr.max(axis=z_ax)
        axes = axes.replace("z", "")

    if keep_channel:
        if "c" not in axes:
            arr = arr[None]
            axes = "c" + axes
        order = [axes.index(a) for a in "cyx"]
    else:
        if "c" in axes:
            c_ax = axes.index("c")
            if arr.shape[c_ax] != 1:
                raise ValueError(f"Expected a single channel, got {arr.shape[c_ax]}")
            arr = arr.take(0, axis=c_ax)
            axes = axes.replace("c", "")
        order = [axes.index(a) for a in "yx"]

    arr = arr.transpose(order)
    return np.asarray(arr.compute() if hasattr(arr, "compute") else arr)


def pixel_size(metadata: dict, level: int) -> tuple[tuple[float, float], str | None]:
    spatial_axes = [ax for ax in metadata["axes"] if ax in "zyx"]
    scales = metadata["scales"][level]
    units = metadata.get("units") or [None] * len(spatial_axes)
    py = float(scales[spatial_axes.index("y")])
    px = float(scales[spatial_axes.index("x")])
    unit = units[spatial_axes.index("x")]
    return (py, px), unit


def match_shape(array: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    """Crop/zero-pad `array` at the high end of each axis to `target_shape`.

    Independently-built multiscale pyramids can round a coarse level's
    dimensions a pixel or two differently even from the same base image.
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


def resolve_channels(specs: list[str], channel_names: list[str]) -> list[int]:
    """Channel indices for `specs`: each an index or an exact channel name."""
    indices = []
    for spec in specs:
        if spec.lstrip("-").isdigit():
            idx = int(spec)
            if not -len(channel_names) <= idx < len(channel_names):
                raise ValueError(f"channel index {idx} out of range {channel_names}")
            indices.append(idx % len(channel_names))
        elif spec in channel_names:
            indices.append(channel_names.index(spec))
        else:
            raise ValueError(f"channel '{spec}' not found in {channel_names}")
    return indices


def column_labels(channel_names: list[str], indices: list[int]) -> list[str]:
    """Channel names usable in column names; a repeated name gets its index appended."""
    labels = [channel_names[i].strip().replace(" ", "_") or f"ch{i}" for i in indices]
    return [f"{lab}_{i}" if labels.count(lab) > 1 else lab for lab, i in zip(labels, indices)]


# --- per-cell morphology + band intensities, and per-image summary --------------


def orientation_degrees(skimage_orientation: float) -> float:
    """Angle of the major axis w.r.t. the image x-axis, in (-90, 90].

    0 = horizontal, positive = counter-clockwise as the image is displayed
    (row 0 at the top). skimage reports the angle to the row axis instead.
    """
    angle = np.degrees(skimage_orientation) - 90.0
    return angle + 180.0 if angle <= -90.0 else angle


def cells_columns(n_regions: int, channel_labels: list[str]) -> list[str]:
    # Fixed so the result always has these columns, even for an image with
    # zero measured cells.
    regions = [f"R{k}" for k in range(1, n_regions + 1)] + ["cell"]
    return (
        [
            "label", "centroid_y_px", "centroid_x_px", "area_px", "area", "perimeter",
            "major_axis_length", "minor_axis_length", "aspect_ratio", "orientation_deg",
        ]
        + [f"R{k}_area_px" for k in range(1, n_regions + 1)]
        + [f"{region}_{label}_mean" for region in regions for label in channel_labels]
    )


def measure_cells(
    labels: np.ndarray, regions: np.ndarray, image: np.ndarray, n_regions: int,
    channels: list[int], channel_labels: list[str], pixel_size_yx: tuple[float, float],
) -> pd.DataFrame:
    """One row per cell that create_endomorph_regions.py painted a band pixel for (a
    cell with no painted band pixel was excluded there for touching the image border):
    morphology, plus the mean of each channel in `channels` in each band R1..RN
    (nested, from `regions`) and over the whole cell mask."""
    py, px = pixel_size_yx
    pixel_area = py * px
    rows = []

    for prop in regionprops(labels):
        r0, c0, r1, c1 = prop.bbox
        cell_mask = labels[r0:r1, c0:c1] == prop.label
        region_crop = regions[r0:r1, c0:c1]
        in_band = (region_crop >= 1) & (region_crop <= n_regions)
        if not (cell_mask & in_band).any():
            continue  # excluded for touching the image border: no band was painted

        crop = image[:, r0:r1, c0:c1]
        masks = {f"R{k}": cell_mask & in_band & (region_crop <= k) for k in range(1, n_regions + 1)}
        masks["cell"] = cell_mask

        minor = prop.axis_minor_length
        row = {
            "label": prop.label,
            "centroid_y_px": prop.centroid[0],
            "centroid_x_px": prop.centroid[1],
            "area_px": int(prop.area),
            "area": prop.area * pixel_area,
            # Perimeter/axes are in pixels; convert assuming isotropic XY sampling.
            "perimeter": prop.perimeter * px,
            "major_axis_length": prop.axis_major_length * px,
            "minor_axis_length": minor * px,
            "aspect_ratio": prop.axis_major_length / minor if minor > 0 else np.nan,
            "orientation_deg": orientation_degrees(prop.orientation),
        }
        for k in range(1, n_regions + 1):
            row[f"R{k}_area_px"] = int(masks[f"R{k}"].sum())

        for region, region_mask in masks.items():
            n = int(region_mask.sum())
            for channel, label in zip(channels, channel_labels):
                values = crop[channel][region_mask]
                row[f"{region}_{label}_mean"] = float(values.mean()) if n else np.nan

        rows.append(row)

    return pd.DataFrame(rows, columns=cells_columns(n_regions, channel_labels))


def summarize_image(labels: np.ndarray, regions: np.ndarray, n_regions: int, n_cells_analysed: int,
                    pixel_size_yx: tuple[float, float]) -> dict:
    """Cell counts and areas. Empty = label N+1 and not imaged = label N+2 of `regions`;
    percentages are of the imaged area."""
    pixel_area = pixel_size_yx[0] * pixel_size_yx[1]
    empty = regions == n_regions + 1
    imaged = regions != n_regions + 2
    labelled = labels > 0
    n_imaged = int(imaged.sum())
    n_empty = int(empty.sum())
    unlabelled = int((~labelled & imaged).sum())
    empty_in_cells = int((empty & labelled).sum())
    pct = (lambda n: 100.0 * n / n_imaged) if n_imaged else (lambda n: np.nan)
    return {
        "n_cells": len(np.unique(labels[labelled])),
        "n_cells_analysed": n_cells_analysed,
        "image_area": labels.size * pixel_area,
        "imaged_area": n_imaged * pixel_area,
        "empty_area": n_empty * pixel_area,
        "empty_area_pct": pct(n_empty),
        "unlabelled_area_pct": pct(unlabelled),
        "undetected_cell_area_pct": pct(unlabelled - (n_empty - empty_in_cells)),
        "empty_in_cells_area_pct": pct(empty_in_cells),
    }


def main():
    args = parse_args()
    n_regions = len(parse_token_list(args.ring_depths))
    if n_regions < 1:
        raise ValueError("--ring-depths needs at least one value")

    intensity_source = mm.ZarrManager(str(args.intensity))
    intensity_md = intensity_source.metadata
    intensity_axes = intensity_md["axes"].lower()
    if args.intensity_resolution_level < 0 or args.intensity_resolution_level >= len(intensity_source.data):
        raise ValueError(
            f"--intensity-resolution-level {args.intensity_resolution_level} is out of "
            f"range for input with {len(intensity_source.data)} level(s)"
        )
    image = extract_plane(
        intensity_source.data[args.intensity_resolution_level], intensity_axes, args.timepoint,
        keep_channel=True,
    )
    pixel_size_yx, unit = pixel_size(intensity_md, args.intensity_resolution_level)
    channel_names = list(intensity_md.get("channel_names") or [f"Channel {i}" for i in range(image.shape[0])])

    labels_source = mm.ZarrManager(str(args.labels))
    labels_md = labels_source.metadata
    labels_axes = labels_md["axes"].lower()
    if args.labels_resolution_level < 0 or args.labels_resolution_level >= len(labels_source.data):
        raise ValueError(
            f"--labels-resolution-level {args.labels_resolution_level} is out of "
            f"range for input with {len(labels_source.data)} level(s)"
        )
    labels_pixel_size_yx, _ = pixel_size(labels_md, args.labels_resolution_level)
    if not np.allclose(pixel_size_yx, labels_pixel_size_yx):
        print(f"WARNING: pixel size differs between --intensity {pixel_size_yx} and "
              f"--labels {labels_pixel_size_yx}; using --intensity's")
    if abs(pixel_size_yx[0] - pixel_size_yx[1]) > 1e-6 * pixel_size_yx[1]:
        print(f"WARNING: anisotropic XY pixels {pixel_size_yx}: perimeter/axis lengths use the X size")

    labels = extract_plane(
        labels_source.data[args.labels_resolution_level], labels_axes, args.timepoint,
        keep_channel=False,
    )
    if not np.issubdtype(labels.dtype, np.integer):
        labels = labels.astype(np.int64)
    if labels.shape != image.shape[1:]:
        print(f"Adjusting labels shape from {labels.shape} to {image.shape[1:]} (pyramid rounding mismatch)")
        labels = match_shape(labels, image.shape[1:])

    regions_source = mm.ZarrManager(str(args.regions))
    regions_axes = regions_source.metadata["axes"].lower()
    regions = extract_plane(regions_source.data[0], regions_axes, 0, keep_channel=False)
    if regions.shape != image.shape[1:]:
        print(f"Adjusting regions shape from {regions.shape} to {image.shape[1:]} (pyramid rounding mismatch)")
        regions = match_shape(regions, image.shape[1:])
    if regions.max() > n_regions + 2:
        raise ValueError(
            f"--regions has label {int(regions.max())} but --ring-depths gives {n_regions} band(s) "
            f"(labels up to {n_regions + 2}): pass the ring depths the regions were made with"
        )

    specs = parse_token_list(args.measure_channels)
    channels = resolve_channels(specs, channel_names) if specs else list(range(image.shape[0]))
    channel_labels = column_labels(channel_names, channels)
    print(f"Measured channels: {dict(zip(channel_labels, channels))}")

    cells = measure_cells(labels, regions, image, n_regions, channels, channel_labels, pixel_size_yx)
    if args.dataset_id:
        cells.insert(0, "dataset_id", args.dataset_id)
    cells.to_csv(args.output, index=False)

    summary = summarize_image(labels, regions, n_regions, len(cells), pixel_size_yx)
    for col in ("area", "aspect_ratio", "perimeter"):
        summary[f"mean_{col}"] = cells[col].mean() if col in cells.columns and len(cells) else np.nan

    py, px = pixel_size_yx
    summary = {
        **({"dataset_id": args.dataset_id} if args.dataset_id else {}),
        **summary,
        "pixel_size_y": py,
        "pixel_size_x": px,
        "unit": unit,
        "n_regions": n_regions,
        "measured_channels": json.dumps(dict(zip(channel_labels, channels))),
    }
    pd.DataFrame([summary]).to_csv(args.output_summary, index=False)

    intensity_source.close()
    labels_source.close()
    regions_source.close()

    print(f"Empty area: {summary['empty_area_pct']:.2f}% of the imaged area")
    print(f"Wrote {len(cells)} cell row(s) to {args.output}")
    print(f"Wrote image summary to {args.output_summary}")


if __name__ == "__main__":
    main()
