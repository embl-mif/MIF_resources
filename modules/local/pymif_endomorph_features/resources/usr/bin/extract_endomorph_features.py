#!/usr/bin/env python3

"""
Measure per-cell morphology and V-cadherin/myosin/actin band intensities
against the R1/R2/R3 edge-distance bands written by
create_endomorph_regions.py, and estimate the cell-free ("empty") area of the
image from the fluorescence signal alone (never from the labels, since
Cellpose can both miss cells and over-reach into cell-free pixels).

A cell is measured if create_endomorph_regions.py painted any band pixel for
it; a cell with no painted pixel was excluded there for touching the image
border.

Writes one CSV of per-cell measurements (morphology + band intensities) and
one CSV of per-image summary metrics (cell counts and empty-area
statistics).

Stripped down from endo_morph's channels.py/empty_area.py/features.py/
run_endo_morph.py (Z:\\people\\gritti\\projects\\Candelas_Adrian\\endo_morph).
"""

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import pymif.microscope_manager as mm
from scipy import ndimage as ndi
from skimage.measure import regionprops

EDGE_THRESHOLD = 0.5  # cell signal is normalized to 0 = background, 1 = cell level

# Keywords are matched case-insensitively against the intensity channel_names.
DEFAULT_KEYWORDS: dict[str, tuple[str, ...]] = {
    "vcad": ("488", "cad"),
    "myosin": ("561", "568", "myo"),
    "actin": ("647", "actin"),
    "nuclei": ("hoechst", "dapi", "405"),
}

# Which markers are averaged in which region.
DEFAULT_MEASUREMENTS: dict[str, tuple[str, ...]] = {
    "R1": ("vcad",),
    "R2": ("myosin", "actin"),
    "R3": ("myosin", "actin"),
    "cell": ("vcad", "myosin", "actin"),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Measure per-cell morphology and band intensities against the "
            "R1/R2/R3 bands of --regions, and the fluorescence-signal empty "
            "area of an image, and write both as CSVs."
        )
    )
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--labels-resolution-level", type=int, default=0)
    parser.add_argument("--timepoint", type=int, default=0)

    ch = parser.add_argument_group("channels (index or exact name; default: keyword match)")
    for marker, keywords in DEFAULT_KEYWORDS.items():
        ch.add_argument(f"--{marker}-channel", default=None,
                        help=f"Keywords used if omitted: {', '.join(keywords)}")

    em = parser.add_argument_group("empty area")
    em.add_argument("--empty-channels", type=str, default="vcad myosin actin",
                    help="Markers whose signal means cells are present, separated by "
                         "comma/semicolon/whitespace. The nuclei stain is excluded by default.")
    em.add_argument("--empty-smoothing", type=float, default=2.0, metavar="UM",
                    help="Gaussian smoothing of the signal before thresholding, in physical units.")
    em.add_argument("--empty-min-area", type=float, default=25.0, metavar="UM2",
                    help="Empty regions smaller than this (physical units^2) are ignored.")
    em.add_argument("--empty-seed-percentile", type=float, default=1.0, metavar="PCT",
                    help="An empty region starts where a pixel is darker than this percentile "
                         "of the cell pixels in every channel at once.")
    em.add_argument("--empty-min-contrast", type=float, default=0.25, metavar="FRACTION",
                    help="A channel contributes to the cell signal only if "
                         "(cell - background) / cell reaches this.")

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


# --- channel resolution --------------------------------------------------------


def resolve_channel(marker: str, spec: str | None, channel_names: list[str]) -> int:
    """Return the channel index for `marker`; `spec` is an index, an exact channel
    name, or None to fall back to keyword matching on the channel names."""
    if spec is not None:
        if spec.lstrip("-").isdigit():
            idx = int(spec)
            if not -len(channel_names) <= idx < len(channel_names):
                raise ValueError(f"{marker}: channel index {idx} out of range {channel_names}")
            return idx % len(channel_names)
        if spec in channel_names:
            return channel_names.index(spec)
        raise ValueError(f"{marker}: channel '{spec}' not found in {channel_names}")

    keywords = DEFAULT_KEYWORDS[marker]
    hits = [
        i for i, name in enumerate(channel_names)
        if any(k in name.lower().replace(" ", "").replace("-", "") for k in keywords)
    ]
    if len(hits) != 1:
        raise ValueError(
            f"{marker}: could not uniquely match keywords {keywords} in channel names "
            f"{channel_names} (matches: {hits}). Pass --{marker}-channel explicitly."
        )
    return hits[0]


def resolve_channels(channel_names: list[str], overrides: dict[str, str | None]) -> dict[str, int]:
    mapping = {m: resolve_channel(m, overrides.get(m), channel_names) for m in DEFAULT_KEYWORDS}
    if len(set(mapping.values())) != len(mapping):
        raise ValueError(f"Several markers map to the same channel: {mapping}")
    return mapping


# --- empty area ----------------------------------------------------------------


@dataclass
class ChannelLevels:
    channel: int
    cell: float
    background: float = np.nan
    used: bool = True

    @property
    def contrast(self) -> float:
        if not self.cell > 0 or not np.isfinite(self.background):
            return 0.0
        return (self.cell - self.background) / self.cell

    @property
    def background_ratio(self) -> float:
        return self.background / self.cell if self.cell > 0 else np.inf


@dataclass
class EmptyArea:
    signal: np.ndarray
    imaged: np.ndarray
    empty: np.ndarray
    seeds: np.ndarray
    threshold: float
    levels: list[ChannelLevels] = field(default_factory=list)

    @property
    def used(self) -> list[int]:
        return [lv.channel for lv in self.levels if lv.used]

    @property
    def background_ratio(self) -> float:
        ratios = [lv.background_ratio for lv in self.levels if lv.used]
        return min(ratios) if ratios else np.inf


def remove_small_regions(mask: np.ndarray, min_px: int) -> np.ndarray:
    regions, n = ndi.label(mask)
    if n == 0:
        return mask
    keep = np.bincount(regions.ravel()) >= min_px
    keep[0] = False
    return keep[regions]


def smooth(values: np.ndarray, weight: np.ndarray, sigma_px: float) -> np.ndarray:
    """Gaussian smoothing that ignores pixels outside `weight` (normalized convolution)."""
    num = ndi.gaussian_filter(values * weight, sigma_px)
    den = ndi.gaussian_filter(weight, sigma_px)
    return num / np.maximum(den, 1e-6)


def grow_from_seeds(signal: np.ndarray, imaged: np.ndarray, seeds: np.ndarray,
                    edge_threshold: float) -> np.ndarray:
    regions, n = ndi.label((signal < edge_threshold) & imaged)
    if n == 0:
        return np.zeros_like(imaged)
    keep = np.zeros(n + 1, bool)
    keep[np.unique(regions[seeds])] = True
    keep[0] = False
    return keep[regions]


def estimate_empty_area(
    image: np.ndarray, labels: np.ndarray, channels: list[int], sigma_px: float,
    min_area_px: float, seed_percentile: float, min_contrast: float,
) -> EmptyArea:
    """See module docstring. `image` is (C, Y, X), `labels` (Y, X)."""
    min_area_px = max(int(round(min_area_px)), 1)

    imaged = ~remove_small_regions(~np.any(image != 0, axis=0), min_area_px)
    weight = imaged.astype(np.float32)

    labelled = labels > 0
    interior = (ndi.distance_transform_edt(labelled) > 2 * sigma_px) & imaged
    if not interior.any():
        interior = labelled & imaged
    reference = interior if interior.any() else imaged

    smoothed = {c: smooth(image[c].astype(np.float32), weight, sigma_px) for c in channels}
    levels = [ChannelLevels(c, float(np.median(smoothed[c][reference]))) for c in channels]

    seeds = imaged.copy()
    for lv in levels:
        seeds &= smoothed[lv.channel] < np.percentile(smoothed[lv.channel][reference], seed_percentile)
    seeds = remove_small_regions(seeds, min_area_px)

    signal = np.zeros(image.shape[1:], np.float32)
    result = EmptyArea(signal, imaged, np.zeros_like(imaged), seeds, 0.0, levels)
    if not seeds.any():
        return result

    for lv in levels:
        lv.background = float(np.median(smoothed[lv.channel][seeds]))
        lv.used = lv.contrast >= min_contrast
        if not lv.used:
            continue
        normalized = np.clip(
            (image[lv.channel] - lv.background) / (lv.cell - lv.background), 0, 1
        ).astype(np.float32)
        np.maximum(signal, smooth(normalized, weight, sigma_px), out=signal)

    if result.used:
        result.threshold = EDGE_THRESHOLD
        result.empty = remove_small_regions(
            grow_from_seeds(signal, imaged, seeds, EDGE_THRESHOLD), min_area_px)
    return result


# --- per-cell morphology + band intensities, and per-image summary --------------


def orientation_degrees(skimage_orientation: float) -> float:
    """Angle of the major axis w.r.t. the image x-axis, in (-90, 90].

    0 = horizontal, positive = counter-clockwise as the image is displayed
    (row 0 at the top). skimage reports the angle to the row axis instead.
    """
    angle = np.degrees(skimage_orientation) - 90.0
    return angle + 180.0 if angle <= -90.0 else angle


def cells_columns(measurements: dict[str, tuple[str, ...]]) -> list[str]:
    # Fixed so the result always has these columns, even for an image with
    # zero measured cells.
    return [
        "label", "centroid_y_px", "centroid_x_px", "area_px", "area", "perimeter",
        "major_axis_length", "minor_axis_length", "aspect_ratio", "orientation_deg",
        "R1_area_px", "R2_area_px", "R3_area_px",
    ] + [f"{region}_{marker}_mean" for region, markers in measurements.items() for marker in markers]


def measure_cells(
    labels: np.ndarray, regions: np.ndarray, image: np.ndarray, channel_map: dict[str, int],
    pixel_size_yx: tuple[float, float], measurements: dict[str, tuple[str, ...]] = DEFAULT_MEASUREMENTS,
) -> pd.DataFrame:
    """One row per cell that create_endomorph_regions.py painted a band pixel for (a
    cell with no painted pixel was excluded there for touching the image border):
    morphology, plus the mean of each requested marker in each region (R1/R2/R3
    bands from `regions`, plus the whole cell mask)."""
    py, px = pixel_size_yx
    pixel_area = py * px
    rows = []

    for prop in regionprops(labels):
        r0, c0, r1, c1 = prop.bbox
        cell_mask = labels[r0:r1, c0:c1] == prop.label
        region_crop = regions[r0:r1, c0:c1]
        if not region_crop[cell_mask].any():
            continue  # excluded for touching the image border: no band was painted

        crop = image[:, r0:r1, c0:c1]
        band_masks = {
            "R1": cell_mask & (region_crop == 1),
            "R2": cell_mask & (region_crop >= 1) & (region_crop <= 2),
            "R3": cell_mask & (region_crop >= 1) & (region_crop <= 3),
            "cell": cell_mask,
        }

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
        for name in ("R1", "R2", "R3"):
            row[f"{name}_area_px"] = int(band_masks[name].sum())

        for region, markers in measurements.items():
            region_mask = band_masks[region]
            n = int(region_mask.sum())
            for marker in markers:
                values = crop[channel_map[marker]][region_mask]
                row[f"{region}_{marker}_mean"] = float(values.mean()) if n else np.nan

        rows.append(row)

    return pd.DataFrame(rows, columns=cells_columns(measurements))


def summarize_image(labels: np.ndarray, n_cells_analysed: int, pixel_size_yx: tuple[float, float],
                    empty: EmptyArea) -> dict:
    pixel_area = pixel_size_yx[0] * pixel_size_yx[1]
    imaged = int(empty.imaged.sum())
    n_empty = int(empty.empty.sum())
    labelled = labels > 0
    unlabelled = int((~labelled & empty.imaged).sum())
    empty_in_cells = int((empty.empty & labelled).sum())
    pct = (lambda n: 100.0 * n / imaged) if imaged else (lambda n: np.nan)
    return {
        "n_cells": len(np.unique(labels[labels > 0])),
        "n_cells_analysed": n_cells_analysed,
        "image_area": labels.size * pixel_area,
        "imaged_area": imaged * pixel_area,
        "empty_area": n_empty * pixel_area,
        "empty_area_pct": pct(n_empty),
        "unlabelled_area_pct": pct(unlabelled),
        "undetected_cell_area_pct": pct(unlabelled - (n_empty - empty_in_cells)),
        "empty_in_cells_area_pct": pct(empty_in_cells),
        "empty_threshold": empty.threshold,
        "empty_background_ratio": empty.background_ratio,
        "empty_seed_area_pct": pct(int(empty.seeds.sum())),
        "empty_channel_levels": json.dumps({
            lv.channel: {"background": round(lv.background, 2), "cell": round(lv.cell, 2), "used": lv.used}
            for lv in empty.levels
        }),
        "empty_channels_used": " ".join(str(c) for c in empty.used),
    }


def main():
    args = parse_args()

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

    overrides = {m: getattr(args, f"{m}_channel") for m in DEFAULT_KEYWORDS}
    channel_map = resolve_channels(channel_names, overrides)
    print(f"Channels: {{{', '.join(f'{m}: {channel_names[i]!r}' for m, i in channel_map.items())}}}")

    cells = measure_cells(labels, regions, image, channel_map, pixel_size_yx)
    if args.dataset_id:
        cells.insert(0, "dataset_id", args.dataset_id)
    cells.to_csv(args.output, index=False)

    empty_channels = parse_token_list(args.empty_channels)
    unknown = [m for m in empty_channels if m not in DEFAULT_KEYWORDS]
    if unknown:
        raise ValueError(f"--empty-channels: unknown marker(s) {unknown}, expected any of {list(DEFAULT_KEYWORDS)}")
    py, px = pixel_size_yx
    empty = estimate_empty_area(
        image, labels, [channel_map[m] for m in empty_channels],
        sigma_px=args.empty_smoothing / px, min_area_px=args.empty_min_area / (py * px),
        seed_percentile=args.empty_seed_percentile, min_contrast=args.empty_min_contrast,
    )
    summary = summarize_image(labels, len(cells), pixel_size_yx, empty)
    for col in ("area", "aspect_ratio", "perimeter"):
        summary[f"mean_{col}"] = cells[col].mean() if col in cells.columns and len(cells) else np.nan

    summary = {
        **({"dataset_id": args.dataset_id} if args.dataset_id else {}),
        **summary,
        "pixel_size_y": py,
        "pixel_size_x": px,
        "unit": unit,
        "empty_channels": " ".join(empty_channels),
        "empty_smoothing_um": args.empty_smoothing,
        "empty_min_area_um2": args.empty_min_area,
        "empty_seed_percentile": args.empty_seed_percentile,
        "empty_min_contrast": args.empty_min_contrast,
        "channels": json.dumps({m: channel_names[i] for m, i in channel_map.items()}),
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
