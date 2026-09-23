#!/usr/bin/env python3

"""
Load a PyMIF label OME-Zarr (e.g. the final Cellpose labels of
conv_mip_cpsam_filt_mask_cpsam) and, for every labelled cell that does not
touch the image border, compute three nested edge-distance bands (R1/R2/R3,
measured inward from the cell edge).

Writes a single-level OME-Zarr label image holding the narrowest band
containing each pixel (1=R1, 2=R2 outside R1, 3=R3 outside R2; 0 elsewhere,
including cells excluded for touching the border -- a cell with no painted
pixel was excluded). All per-cell morphology and marker intensities are
computed downstream, from this zarr, by extract_endomorph_features.py.

Stripped down from endo_morph's regions.py/zarr_io.py
(Z:\\people\\gritti\\projects\\Candelas_Adrian\\endo_morph).
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pymif.microscope_manager as mm
from scipy import ndimage as ndi
from skimage.measure import regionprops

REGION_NAMES = ("R1", "R2", "R3")

# NGFF requires UCUM-style unit names; pymif reports e.g. "µm".
NGFF_UNITS = {
    "µm": "micrometer", "μm": "micrometer", "um": "micrometer", "micron": "micrometer",
    "nm": "nanometer", "mm": "millimeter",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute nested R1/R2/R3 edge-distance bands for every "
            "non-border-touching cell of a label OME-Zarr, and write them "
            "as a label OME-Zarr."
        )
    )
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--resolution-level", type=int, default=0)
    parser.add_argument("--timepoint", type=int, default=0)
    parser.add_argument(
        "--ring-depths",
        type=str,
        default="6 9 24",
        help=(
            "Three band depths (R1 R2 R3), from the cell edge in pixels of "
            "--resolution-level, separated by comma/semicolon/whitespace. "
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
    return parser.parse_args()


def parse_ring_depths(text: str) -> tuple[float, float, float]:
    tokens = [token for token in re.split(r"[;,\s]+", text.strip()) if token]
    if len(tokens) != 3:
        raise ValueError(f"--ring-depths must have exactly 3 values (R1 R2 R3), got {text!r}")
    return tuple(float(token) for token in tokens)


def extract_plane(level_data, axes: str, t: int) -> np.ndarray:
    """Select timepoint `t`, max-project Z if present, drop a size-1 channel axis, and
    return a (Y, X) array."""
    axes = axes.lower()
    indexer = [t if ax == "t" else slice(None) for ax in axes]
    arr = level_data[tuple(indexer)]
    axes = axes.replace("t", "")

    if "z" in axes:
        z_ax = axes.index("z")
        arr = arr.max(axis=z_ax)
        axes = axes.replace("z", "")

    if "c" in axes:
        c_ax = axes.index("c")
        if arr.shape[c_ax] != 1:
            raise ValueError(f"--labels has {arr.shape[c_ax]} channel(s), expected a single label channel")
        arr = arr.take(0, axis=c_ax)
        axes = axes.replace("c", "")

    order = [axes.index(a) for a in "yx"]
    arr = arr.transpose(order)
    return np.asarray(arr.compute() if hasattr(arr, "compute") else arr)


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


def main():
    args = parse_args()
    depths = parse_ring_depths(args.ring_depths)

    source = mm.ZarrManager(str(args.labels))
    metadata = source.metadata
    axes = metadata["axes"].lower()
    if args.resolution_level < 0 or args.resolution_level >= len(source.data):
        raise ValueError(
            f"--resolution-level {args.resolution_level} is out of range for "
            f"input with {len(source.data)} level(s)"
        )

    labels = extract_plane(source.data[args.resolution_level], axes, args.timepoint)
    if not np.issubdtype(labels.dtype, np.integer):
        labels = labels.astype(np.int64)

    spatial_axes = [ax for ax in axes if ax in "zyx"]
    scales = metadata["scales"][args.resolution_level]
    units = metadata.get("units") or [None] * len(spatial_axes)
    py = float(scales[spatial_axes.index("y")])
    px = float(scales[spatial_axes.index("x")])
    unit = units[spatial_axes.index("x")]

    h, w = labels.shape
    region_image = np.zeros(labels.shape, np.uint8)
    n_total = 0
    n_measured = 0

    for prop in regionprops(labels):
        n_total += 1
        r0, c0, r1, c1 = prop.bbox
        # bbox is [r0, r1) x [c0, c1); margin 0 means touching the first/last pixel.
        if min(r0, c0, h - r1, w - c1) <= args.border_margin:
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

    print(
        f"{n_measured}/{n_total} cell(s) measured (excluded {n_total - n_measured} "
        f"touching the border within {args.border_margin}px)"
    )
    print(f"Bands (px from edge): {dict(zip(REGION_NAMES, depths))}")
    print(f"Wrote region bands to {args.output}")


if __name__ == "__main__":
    main()
