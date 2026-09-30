#!/usr/bin/env python3

"""
Load a PyMIF label OME-Zarr (e.g. the final Cellpose labels of
conv_mip_cpsam_filt_mask_cpsam) and its MIP intensity OME-Zarr and write ONE
label OME-Zarr describing where things are:

  1=R1, 2=R2 outside R1, 3=R3 outside R2   nested edge-distance bands (measured
                                           inward from the cell edge) of every
                                           labelled cell that does not touch the
                                           image border
  4=empty                                  cell-free area, estimated from the
                                           fluorescence signal alone (never from
                                           the labels, since Cellpose can both
                                           miss cells and over-reach into
                                           cell-free pixels)
  5=not imaged                             zero in every channel (e.g. padding)
  0=everything else, including cells excluded for touching the border (a cell
    with no painted band pixel was excluded)

Bands keep priority: 4 and 5 are painted only on pixels without a band. Because
that hides empty pixels that fall inside cell masks, all the image-level
empty-area statistics are computed here from the full masks and written to a
JSON file (--output-empty-stats) for extract_endomorph_features.py to add to its
per-image summary.

The zarr is written at the resolution of the intensity image
(--intensity-resolution-level); the labels are cropped/zero-padded to it if
independent pyramids rounded the shape differently.

Stripped down from endo_morph's regions.py/empty_area.py/zarr_io.py
(Z:\\people\\gritti\\projects\\Candelas_Adrian\\endo_morph).
"""

import argparse
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pymif.microscope_manager as mm
from scipy import ndimage as ndi
from skimage.measure import regionprops

REGION_NAMES = ("R1", "R2", "R3")
QC_EMPTY = 4
QC_NOT_IMAGED = 5

EDGE_THRESHOLD = 0.5  # cell signal is normalized to 0 = background, 1 = cell level

# Keywords are matched case-insensitively against the intensity channel_names.
DEFAULT_KEYWORDS: dict[str, tuple[str, ...]] = {
    "vcad": ("488", "cad"),
    "myosin": ("561", "568", "myo"),
    "actin": ("647", "actin"),
    "nuclei": ("hoechst", "dapi", "405"),
}

# NGFF requires UCUM-style unit names; pymif reports e.g. "µm".
NGFF_UNITS = {
    "µm": "micrometer", "μm": "micrometer", "um": "micrometer", "micron": "micrometer",
    "nm": "nanometer", "mm": "millimeter",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute nested R1/R2/R3 edge-distance bands for every "
            "non-border-touching cell of a label OME-Zarr, plus the empty and "
            "not-imaged areas of the intensity image, and write them as one "
            "label OME-Zarr."
        )
    )
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-empty-stats", required=True, type=Path,
                        help="JSON file with the image-level empty-area statistics.")
    parser.add_argument("--resolution-level", type=int, default=0,
                        help="Resolution level of --labels.")
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--timepoint", type=int, default=0)
    parser.add_argument(
        "--ring-depths",
        type=str,
        default="6 9 24",
        help=(
            "Three band depths (R1 R2 R3), from the cell edge in pixels of "
            "the intensity resolution level, separated by comma/semicolon/whitespace. "
            "All bands start at the edge, so they are nested."
        ),
    )
    parser.add_argument(
        "--border-margin",
        type=int,
        default=10,
        help=(
            "Cells with any pixel within this many pixels of the image "
            "border are considered cut by the field of view and excluded."
        ),
    )

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
    return parser.parse_args()


def parse_token_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [token for token in re.split(r"[;,\s]+", text.strip()) if token]


def parse_ring_depths(text: str) -> tuple[float, float, float]:
    tokens = parse_token_list(text)
    if len(tokens) != 3:
        raise ValueError(f"--ring-depths must have exactly 3 values (R1 R2 R3), got {text!r}")
    return tuple(float(token) for token in tokens)


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


# --- bands ---------------------------------------------------------------------


def edge_distance(cell_mask: np.ndarray) -> np.ndarray:
    """Distance (px) of every in-mask pixel to the nearest pixel outside the mask.

    The mask must be surrounded by background, which always holds here: cells
    touching the image border are excluded before this is called.
    """
    dist = ndi.distance_transform_edt(cell_mask)
    dist[~cell_mask] = 0
    return dist


def band_masks(dist: np.ndarray, depths: tuple[float, ...]) -> list[np.ndarray]:
    """One mask per depth, each covering ``0 < dist <= depth`` (nested bands)."""
    return [(dist > 0) & (dist <= d) for d in depths]


def paint_bands(labels: np.ndarray, depths: tuple[float, ...], border_margin: int) -> tuple[np.ndarray, int, int]:
    """Band image (0/1/2/3) plus the number of cells measured and in total."""
    h, w = labels.shape
    region_image = np.zeros(labels.shape, np.uint8)
    n_total = 0
    n_measured = 0

    for prop in regionprops(labels):
        n_total += 1
        r0, c0, r1, c1 = prop.bbox
        # bbox is [r0, r1) x [c0, c1); margin 0 means touching the first/last pixel.
        if min(r0, c0, h - r1, w - c1) <= border_margin:
            continue  # cut by the image border: not measured

        # Crop with a 1-px background margin so the mask edge is detected everywhere.
        sr0, sc0, sr1, sc1 = r0 - 1, c0 - 1, r1 + 1, c1 + 1
        mask = labels[sr0:sr1, sc0:sc1] == prop.label

        dist = edge_distance(mask)
        bands = dict(zip(REGION_NAMES, band_masks(dist, depths)))
        n_measured += 1

        target = region_image[sr0:sr1, sc0:sc1]
        # Paint the widest band first so narrower bands stay visible on top.
        order = sorted(enumerate(REGION_NAMES, start=1), key=lambda x: -depths[x[0] - 1])
        for i, name in order:
            target[bands[name]] = i

    return region_image, n_measured, n_total


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
    image: np.ndarray,
    labels: np.ndarray,
    channels: list[int],
    sigma_px: float = 2.0,
    min_area_px: int = 25,
    seed_percentile: float = 1.0,
    min_contrast: float = 0.25,
) -> EmptyArea:
    """`image` is (C, Y, X), `labels` (Y, X). Empty pixels are dark in *every* marker
    channel, cells are dark in some channel almost everywhere; labels are only used to
    scale the cell level."""
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


def empty_area_stats(labels: np.ndarray, pixel_size_yx: tuple[float, float], empty: EmptyArea) -> dict:
    """Image-level empty-area statistics, from the full masks (before band priority)."""
    pixel_area = pixel_size_yx[0] * pixel_size_yx[1]
    imaged = int(empty.imaged.sum())
    n_empty = int(empty.empty.sum())
    labelled = labels > 0
    unlabelled = int((~labelled & empty.imaged).sum())
    empty_in_cells = int((empty.empty & labelled).sum())
    pct = (lambda n: 100.0 * n / imaged) if imaged else (lambda n: float("nan"))
    return {
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
    depths = parse_ring_depths(args.ring_depths)

    intensity_source = mm.ZarrManager(str(args.intensity))
    intensity_md = intensity_source.metadata
    if args.intensity_resolution_level < 0 or args.intensity_resolution_level >= len(intensity_source.data):
        raise ValueError(
            f"--intensity-resolution-level {args.intensity_resolution_level} is out of "
            f"range for input with {len(intensity_source.data)} level(s)"
        )
    image = extract_plane(
        intensity_source.data[args.intensity_resolution_level], intensity_md["axes"], args.timepoint,
        keep_channel=True,
    )
    (py, px), unit = pixel_size(intensity_md, args.intensity_resolution_level)
    channel_names = list(intensity_md.get("channel_names") or [f"Channel {i}" for i in range(image.shape[0])])

    source = mm.ZarrManager(str(args.labels))
    if args.resolution_level < 0 or args.resolution_level >= len(source.data):
        raise ValueError(
            f"--resolution-level {args.resolution_level} is out of range for "
            f"input with {len(source.data)} level(s)"
        )
    labels = extract_plane(
        source.data[args.resolution_level], source.metadata["axes"], args.timepoint, keep_channel=False,
    )
    if not np.issubdtype(labels.dtype, np.integer):
        labels = labels.astype(np.int64)
    if labels.shape != image.shape[1:]:
        print(f"Adjusting labels shape from {labels.shape} to {image.shape[1:]} (pyramid rounding mismatch)")
        labels = match_shape(labels, image.shape[1:])

    region_image, n_measured, n_total = paint_bands(labels, depths, args.border_margin)

    overrides = {m: getattr(args, f"{m}_channel") for m in DEFAULT_KEYWORDS}
    channel_map = resolve_channels(channel_names, overrides)
    print(f"Channels: {{{', '.join(f'{m}: {channel_names[i]!r}' for m, i in channel_map.items())}}}")
    empty_channels = parse_token_list(args.empty_channels)
    unknown = [m for m in empty_channels if m not in DEFAULT_KEYWORDS]
    if unknown:
        raise ValueError(f"--empty-channels: unknown marker(s) {unknown}, expected any of {list(DEFAULT_KEYWORDS)}")

    empty = estimate_empty_area(
        image, labels, [channel_map[m] for m in empty_channels],
        sigma_px=args.empty_smoothing / px, min_area_px=args.empty_min_area / (py * px),
        seed_percentile=args.empty_seed_percentile, min_contrast=args.empty_min_contrast,
    )
    stats = {
        **empty_area_stats(labels, (py, px), empty),
        "empty_channels": " ".join(empty_channels),
        "empty_smoothing_um": args.empty_smoothing,
        "empty_min_area_um2": args.empty_min_area,
        "empty_seed_percentile": args.empty_seed_percentile,
        "empty_min_contrast": args.empty_min_contrast,
    }
    args.output_empty_stats.write_text(json.dumps(stats, indent=2))

    # Bands keep priority: empty / not imaged only fill pixels without a band.
    region_image[empty.empty & (region_image == 0)] = QC_EMPTY
    region_image[~empty.imaged & (region_image == 0)] = QC_NOT_IMAGED

    out_unit = NGFF_UNITS.get(unit, unit) or "micrometer"
    out_metadata = {
        "axes": "yx",
        "scales": [(py, px)],
        "units": (out_unit, out_unit),
        "data_type": "label",
        "name": Path(str(args.output)).name,
    }
    mm.ArrayManager(region_image, metadata=out_metadata, chunks=(1024, 1024)).to_zarr(
        str(args.output), ngff_version="0.5", overwrite=True
    )

    source.close()
    intensity_source.close()

    print(
        f"{n_measured}/{n_total} cell(s) measured (excluded {n_total - n_measured} "
        f"touching the border within {args.border_margin}px)"
    )
    print(f"Bands (px from edge): {dict(zip(REGION_NAMES, depths))}")
    print(f"Empty area: {stats['empty_area_pct']:.2f}% of the imaged area")
    print(f"Wrote region bands + empty area to {args.output}")


if __name__ == "__main__":
    main()
