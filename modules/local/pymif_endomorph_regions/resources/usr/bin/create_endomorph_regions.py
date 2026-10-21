#!/usr/bin/env python3

"""
Load a PyMIF label OME-Zarr (e.g. the final Cellpose labels) and its intensity
OME-Zarr (e.g. the MIP) and write ONE label OME-Zarr describing where things are.
With N ring depths (any N >= 1):

  1..N     nested edge-distance bands R1..RN, measured inward from the cell edge,
           of every labelled cell that does not touch the image border. A pixel
           holds the narrowest band containing it (1=R1, 2=R2 outside R1, ...).
  N+1      empty (gaps): cell-free area, dark in ALL --empty-channels (the
           cell-body markers, e.g. VE-cadherin and actin) and not a host nucleus
           (--empty-nucleus-channel of --empty-nucleus-zarr, e.g.
           HOECHST_33342_clean of the clean-channel zarr). Decided by the
           fluorescence alone, never by the labels, since Cellpose can both miss
           cells and over-reach into cell-free pixels
  N+2      not imaged: zero in every channel (e.g. padding)
  0        everything else, including cells excluded for touching the border (a
           cell with no painted band pixel was excluded)

Empty and not-imaged pixels overwrite the bands. All measurements (per-cell
features and the per-image summary, including the empty area) are computed
downstream, from this zarr, by extract_endomorph_features.py. A JSON with QC
diagnostics of the empty-area estimation (channel levels, threshold, ...) is
written to --output-empty-qc.

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

NOT_IMAGED_MIN_AREA_UM2 = 25.0  # all-zero regions at least this large count as not imaged

# NGFF requires UCUM-style unit names; pymif reports e.g. "µm".
NGFF_UNITS = {
    "µm": "micrometer", "μm": "micrometer", "um": "micrometer", "micron": "micrometer",
    "nm": "nanometer", "mm": "millimeter",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute nested R1..RN edge-distance bands for every "
            "non-border-touching cell of a label OME-Zarr, plus the empty "
            "(N+1) and not-imaged (N+2) areas of the intensity image, and "
            "write them as one label OME-Zarr."
        )
    )
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--intensity", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-empty-qc", required=True, type=Path,
                        help="JSON file with QC diagnostics of the empty-area estimation.")
    parser.add_argument("--resolution-level", type=int, default=0,
                        help="Resolution level of --labels.")
    parser.add_argument("--intensity-resolution-level", type=int, default=0)
    parser.add_argument("--timepoint", type=int, default=0)
    parser.add_argument(
        "--ring-depths",
        type=str,
        default="6 9 24",
        help=(
            "Band depths (R1 R2 ... RN, any number >= 1, strictly increasing), from the "
            "cell edge in pixels of the intensity resolution level, separated by "
            "comma/semicolon/whitespace. All bands start at the edge, so they are nested."
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

    em = parser.add_argument_group("empty area")
    em.add_argument("--empty-channels", default=None, type=str,
                    help="Intensity channels (index or exact name, separated by "
                         "comma/semicolon/whitespace) whose signal means cells are present. "
                         "A gap must be dark in ALL of them. Default: all channels except the nuclear "
                         "one. Use the cell-body markers only (e.g. VE-cadherin + actin); leave out "
                         "channels that only stain part of the cell or cannot tell a gap from cytoplasm.")
    em.add_argument("--empty-nucleus-zarr", type=Path, default=None,
                    help="OME-Zarr holding the nuclear channel, e.g. the clean-channel zarr. Default: --intensity.")
    em.add_argument("--empty-nucleus-channel", default="auto", type=str,
                    help="Nuclear-stain channel of --empty-nucleus-zarr (index or exact name, e.g. "
                         "HOECHST_33342_clean). Its nuclei are never gap. It, or its source channel "
                         "(the name without '_clean'), is removed from --empty-channels. 'auto' (default): "
                         "the channel whose name contains 'hoechst' or 'dapi' (case-insensitive), "
                         "preferring '..._clean'; 'none': no nucleus veto.")
    em.add_argument("--empty-nucleus-factor", type=float, default=8.0, metavar="X",
                    help="A pixel is nucleus when its smoothed nuclear signal exceeds this multiple "
                         "of the image median of that channel. A background-free '_clean' channel "
                         "(median 0) uses its non-zero pixels instead.")
    em.add_argument("--empty-smoothing", type=float, default=0.3, metavar="UM",
                    help="Gaussian smoothing of the signal before thresholding, in physical units. "
                         "Keep it small (0.25-0.5): more smoothing blurs small holes into their bright rim.")
    em.add_argument("--empty-min-area", type=float, default=3.0, metavar="UM2",
                    help="Empty regions smaller than this (physical units^2) are ignored.")
    em.add_argument("--empty-core-min-area", type=float, default=1.5, metavar="UM2",
                    help="A gap must contain a core at background of at least this area (physical units^2).")
    em.add_argument("--empty-seed-percentile", type=float, default=1.0, metavar="PCT",
                    help="An empty region starts where a pixel is darker than this percentile "
                         "of the cell pixels in every channel at once.")
    em.add_argument("--empty-min-contrast", type=float, default=0.25, metavar="FRACTION",
                    help="A channel contributes to the cell signal only if "
                         "(cell - background) / cell reaches this.")
    em.add_argument("--empty-cell-percentile", type=float, default=10.0, metavar="PCT",
                    help="Cell level of a channel = this (low) percentile of the labelled-cell interior. "
                         "A lower value amplifies the signal (more gain): dimmer cytoplasm counts as cell.")
    em.add_argument("--empty-core-threshold", type=float, default=0.2, metavar="0-1",
                    help="A gap must contain a core where the normalized cell signal is below this "
                         "(0 = background, 1 = cell level); the gap extends to --empty-edge-threshold.")
    em.add_argument("--empty-edge-threshold", type=float, default=0.65, metavar="0-1",
                    help="The gap extends to where the normalized cell signal reaches this value.")
    em.add_argument("--empty-max-floor-ratio", type=float, default=0.55, metavar="RATIO",
                    help="No gap is reported unless the gap floor is below this fraction of the darkest "
                         "cell pixels (guards against inventing gaps in a confluent field).")
    return parser.parse_args()


def parse_token_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [token for token in re.split(r"[;,\s]+", text.strip()) if token]


def parse_ring_depths(text: str) -> tuple[float, ...]:
    depths = tuple(float(token) for token in parse_token_list(text))
    if not depths:
        raise ValueError("--ring-depths needs at least one value")
    if any(d <= 0 for d in depths) or any(a >= b for a, b in zip(depths, depths[1:])):
        raise ValueError(f"--ring-depths must be positive and strictly increasing, got {text!r}")
    if len(depths) + 2 > 255:
        raise ValueError("--ring-depths: too many bands for an 8-bit label image")
    return depths


def find_nucleus_channel(spec: str | None, channel_names: list[str]) -> int | None:
    """Index of the nuclear-stain channel: an explicit index/name, 'auto' (the channel named
    like Hoechst/DAPI, preferring '..._clean'), or None ('none'/blank)."""
    tokens = parse_token_list(spec)
    if not tokens or tokens[0].lower() == "none":
        return None
    if tokens[0].lower() == "auto":
        matches = [i for i, name in enumerate(channel_names) if re.search(r"hoechst|dapi", name, re.IGNORECASE)]
        clean = [i for i in matches if "clean" in channel_names[i].lower()]
        return (clean or matches or [None])[0]
    return resolve_channels(tokens[:1], channel_names)[0]


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


# --- bands ---------------------------------------------------------------------


def edge_distance(cell_mask: np.ndarray) -> np.ndarray:
    """Distance (px) of every in-mask pixel to the nearest pixel outside the mask.

    The mask must be surrounded by background, which always holds here: cells
    touching the image border are excluded before this is called.
    """
    dist = ndi.distance_transform_edt(cell_mask)
    dist[~cell_mask] = 0
    return dist


def paint_bands(labels: np.ndarray, depths: tuple[float, ...], border_margin: int) -> tuple[np.ndarray, int, int]:
    """Band image (0, 1..N) plus the number of cells measured and in total."""
    h, w = labels.shape
    region_image = np.zeros(labels.shape, np.uint8)
    n_total = 0
    n_measured = 0

    # Paint the widest band first so narrower bands stay visible on top.
    order = sorted(enumerate(depths, start=1), key=lambda x: -x[1])

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
        n_measured += 1

        target = region_image[sr0:sr1, sc0:sc1]
        for i, depth in order:
            target[(dist > 0) & (dist <= depth)] = i

    return region_image, n_measured, n_total


# --- empty area ----------------------------------------------------------------


@dataclass
class ChannelLevels:
    channel: int
    cell: float
    background: float = np.nan
    used: bool = True
    seed_level: float = np.nan  # smoothed value below which a cell pixel is in the darkest seed_percentile

    @property
    def floor_ratio(self) -> float:
        """Background (gap) level relative to the darkest cell pixels; << 1 for a real gap."""
        if not self.seed_level > 0 or not np.isfinite(self.background):
            return np.inf
        return self.background / self.seed_level

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
    nucleus: np.ndarray | None = None  # pixels vetoed as nucleus (never gap)
    gap_found: bool = False

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


def fill_small_holes(mask: np.ndarray, max_hole_px: int) -> np.ndarray:
    """Fill holes of `mask` smaller than `max_hole_px` pixels (specks inside a gap)."""
    holes = ndi.binary_fill_holes(mask) & ~mask
    big = remove_small_regions(holes, max_hole_px + 1)
    return mask | (holes & ~big)


def estimate_empty_area(
    image: np.ndarray,
    labels: np.ndarray,
    channels: list[int],
    sigma_px: float = 1.6,
    min_area_px: int = 87,
    seed_percentile: float = 1.0,
    min_contrast: float = 0.25,
    nucleus_channel: int | None = None,
    nucleus_factor: float = 8.0,
    cell_percentile: float = 10.0,
    core_threshold: float = 0.2,
    edge_threshold: float = 0.65,
    core_min_px: int = 43,
    max_floor_ratio: float = 0.55,
    imaged_min_px: int = 723,
) -> EmptyArea:
    """`image` is (C, Y, X), `labels` (Y, X). Gap = cell-free area, from the fluorescence of the
    cell-body markers `channels` (e.g. VE-cadherin + actin); a gap must be dark in ALL of them.

    1. Nuclei (`nucleus_channel` > `nucleus_factor` x its image median) are never gap: a nucleus is
       naturally dark in VE-cadherin/actin.
    2. Each channel is lightly smoothed (`sigma_px`, small so small holes survive).
    3. Seeds = pixels in the darkest `seed_percentile` of the labelled-cell interior in every channel
       at once; the background (gap floor) of a channel is its median inside the seeds. No seed, or a
       floor that is not clearly darker than the darkest cells (`max_floor_ratio`), means no gap
       (a confluent field).
    4. Gain: the cell level is a LOW percentile (`cell_percentile`) of the cell interior, so ordinary
       cytoplasm saturates at 1 and only background stays low:
       signal = max over channels of clip((smoothed - floor) / (cell - floor), 0, 1); nuclei = 1.
    5. Hysteresis: a gap is a region of signal < `edge_threshold` (the edge; the outline moves
       < 1 px between 0.5 and 0.8, the area < 10 %) containing a core of signal < `core_threshold`
       of at least `core_min_px`, i.e. truly at background, not merely dim.
    Labels are only used to scale the cell level; they never decide what is empty.
    """
    min_area_px = max(int(round(min_area_px)), 1)
    core_min_px = max(int(round(core_min_px)), 1)

    imaged = ~remove_small_regions(~np.any(image != 0, axis=0), max(int(round(imaged_min_px)), 1))
    weight = imaged.astype(np.float32)

    nucleus = np.zeros_like(imaged)
    if nucleus_channel is not None:
        raw_nuc = image[nucleus_channel]
        if np.median(raw_nuc[imaged]) > 0:
            hoechst = smooth(raw_nuc.astype(np.float32), weight, sigma_px)
            nucleus = imaged & (hoechst > nucleus_factor * np.median(hoechst[imaged]))
        else:
            # Background-free "clean" nuclear channel (zero outside the host nuclei): the nuclei are the
            # non-zero pixels. A multiple of the median would be 0 here and smoothing would add a halo.
            nucleus = imaged & (raw_nuc > 0)
    candidates = imaged & ~nucleus

    labelled = labels > 0
    interior = (ndi.distance_transform_edt(labelled) > 2 * sigma_px) & candidates
    if not interior.any():
        interior = labelled & candidates
    reference = interior if interior.any() else candidates

    smoothed = {c: smooth(image[c].astype(np.float32), weight, sigma_px) for c in channels}
    levels = [ChannelLevels(c, float(np.percentile(smoothed[c][reference], cell_percentile)),
                            seed_level=float(np.percentile(smoothed[c][reference], seed_percentile)))
              for c in channels]

    seeds = candidates.copy()
    for lv in levels:
        seeds &= smoothed[lv.channel] < lv.seed_level
    seeds = remove_small_regions(seeds, core_min_px)

    signal = np.zeros(image.shape[1:], np.float32)
    result = EmptyArea(signal, imaged, np.zeros_like(imaged), seeds, 0.0, levels, nucleus)
    if not seeds.any():
        return result

    parts = []
    for lv in levels:
        lv.background = float(np.median(smoothed[lv.channel][seeds]))
        lv.used = lv.contrast >= min_contrast
        if lv.used:
            parts.append(np.clip((smoothed[lv.channel] - lv.background) / (lv.cell - lv.background), 0, 1))
    if not parts or min(lv.floor_ratio for lv in levels if lv.used) > max_floor_ratio:
        return result

    signal[:] = np.maximum.reduce(parts)
    signal[nucleus] = 1.0

    core = remove_small_regions((signal < core_threshold) & imaged, core_min_px)
    result.seeds = core
    result.threshold = edge_threshold
    empty = grow_from_seeds(signal, imaged, core, edge_threshold)
    empty = fill_small_holes(empty, core_min_px)
    result.empty = remove_small_regions(empty, min_area_px)
    result.gap_found = bool(result.empty.any())
    return result


def empty_area_qc(empty: EmptyArea, channel_names: list[str]) -> dict:
    """QC diagnostics of the estimation (what the empty area is not enough to tell)."""
    imaged = int(empty.imaged.sum())
    return {
        "empty_threshold": empty.threshold,
        "empty_background_ratio": empty.background_ratio,
        "empty_seed_area_pct": 100.0 * int(empty.seeds.sum()) / imaged if imaged else float("nan"),
        "empty_channels_used": [channel_names[c] for c in empty.used],
        "empty_channel_levels": {
            channel_names[lv.channel]: {
                "background": round(lv.background, 2), "cell": round(lv.cell, 2), "used": lv.used,
                "floor_ratio": round(lv.floor_ratio, 3) if np.isfinite(lv.floor_ratio) else None,
            }
            for lv in empty.levels
        },
        "empty_gap_found": empty.gap_found,
        "empty_n_gaps": int(ndi.label(empty.empty)[1]),
        "empty_nucleus_veto_area_pct": (
            100.0 * int(empty.nucleus.sum()) / imaged if imaged and empty.nucleus is not None else None
        ),
    }


def main():
    args = parse_args()
    depths = parse_ring_depths(args.ring_depths)
    n_regions = len(depths)

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

    # Nuclear channel: from --empty-nucleus-zarr (e.g. the clean-channel zarr) or from --intensity.
    nuc_image, nuc_names = image, channel_names
    if args.empty_nucleus_zarr:
        nuc_source = mm.ZarrManager(str(args.empty_nucleus_zarr))
        nuc_image = extract_plane(
            nuc_source.data[args.intensity_resolution_level], nuc_source.metadata["axes"], args.timepoint,
            keep_channel=True,
        )
        nuc_names = list(nuc_source.metadata.get("channel_names")
                         or [f"Channel {i}" for i in range(nuc_image.shape[0])])
        nuc_source.close()
    nucleus_idx = find_nucleus_channel(args.empty_nucleus_channel, nuc_names)
    nucleus_name = nuc_names[nucleus_idx] if nucleus_idx is not None else None

    specs = parse_token_list(args.empty_channels)
    empty_channels = resolve_channels(specs, channel_names) if specs else list(range(image.shape[0]))
    nucleus_channel = None
    if nucleus_idx is not None:
        # The nuclear stain (or the source of a "_clean" channel) is never a gap channel.
        source_name = re.sub(r"_clean$", "", nucleus_name)
        empty_channels = [c for c in empty_channels if channel_names[c] not in (nucleus_name, source_name)]
        # Appended after the intensity channels so the channel indices above stay valid.
        image = np.concatenate([image, match_shape(nuc_image[nucleus_idx], image.shape[1:])[None]])
        nucleus_channel = image.shape[0] - 1
        print(f"Empty-area nucleus veto channel: {nucleus_name}")
    if not empty_channels:
        raise ValueError("--empty-channels needs at least one channel other than the nuclear one")
    print(f"Empty-area channels: {[channel_names[c] for c in empty_channels]}")

    empty = estimate_empty_area(
        image, labels, empty_channels,
        sigma_px=args.empty_smoothing / px, min_area_px=args.empty_min_area / (py * px),
        seed_percentile=args.empty_seed_percentile, min_contrast=args.empty_min_contrast,
        nucleus_channel=nucleus_channel, nucleus_factor=args.empty_nucleus_factor,
        cell_percentile=args.empty_cell_percentile, core_threshold=args.empty_core_threshold,
        edge_threshold=args.empty_edge_threshold,
        core_min_px=args.empty_core_min_area / (py * px), max_floor_ratio=args.empty_max_floor_ratio,
        imaged_min_px=NOT_IMAGED_MIN_AREA_UM2 / (py * px),
    )
    qc = {
        **empty_area_qc(empty, channel_names),
        "empty_channels": [channel_names[c] for c in empty_channels],
        "empty_nucleus_channel": nucleus_name,
        "empty_nucleus_factor": args.empty_nucleus_factor,
        "empty_smoothing_um": args.empty_smoothing,
        "empty_min_area_um2": args.empty_min_area,
        "empty_core_min_area_um2": args.empty_core_min_area,
        "empty_seed_percentile": args.empty_seed_percentile,
        "empty_min_contrast": args.empty_min_contrast,
        "empty_cell_percentile": args.empty_cell_percentile,
        "empty_core_threshold": args.empty_core_threshold,
        "empty_edge_threshold": args.empty_edge_threshold,
        "empty_max_floor_ratio": args.empty_max_floor_ratio,
    }
    args.output_empty_qc.write_text(json.dumps(qc, indent=2))

    region_image[empty.empty] = n_regions + 1
    region_image[~empty.imaged] = n_regions + 2

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
    print(f"Bands (px from edge): { {f'R{i}': d for i, d in enumerate(depths, start=1)} }")
    print(f"Labels: 1..{n_regions}=bands, {n_regions + 1}=empty, {n_regions + 2}=not imaged")
    print(f"Wrote region bands + empty area to {args.output}")


if __name__ == "__main__":
    main()
