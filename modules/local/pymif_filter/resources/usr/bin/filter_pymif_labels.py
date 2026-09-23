#!/usr/bin/env python3

"""
Filter an existing PyMIF/OME-Zarr label store against a precomputed
per-label measurements CSV (e.g. from extract_pymif_features.py): keep only
the labels whose row satisfies the requested area and/or other property
limits, zero every other label out of every pyramid level, relabel the kept
objects to a contiguous 1..N id space (consistently across every pyramid
level, and per timepoint if the labels have a time axis), and write both the
filtered/relabeled label OME-Zarr and a CSV holding only the kept rows (with
the new "label" alongside the original id, as "original_label").
"""

import argparse
import copy
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import pymif.microscope_manager as mm


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter a label OME-Zarr using min/max area and/or other property "
            "limits measured in a companion CSV, and write the filtered "
            "labels and CSV."
        )
    )
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--min-area", type=float, default=None)
    parser.add_argument("--max-area", type=float, default=None)
    parser.add_argument(
        "--filters",
        type=str,
        default=None,
        help=(
            "Semicolon-separated 'property:min:max' bounds on any column of "
            "--csv. Either bound may be left empty, e.g. "
            "'solidity:0.8:;mean_intensity_ch1:100:'."
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


def filter_rows(
    df: pd.DataFrame,
    min_area: float | None,
    max_area: float | None,
    filters: dict[str, tuple[float | None, float | None]],
) -> pd.DataFrame:
    keep = pd.Series(True, index=df.index)

    if min_area is not None or max_area is not None:
        if "area" not in df.columns:
            raise ValueError(
                "--min-area/--max-area was given but --csv has no 'area' "
                f"column. Available columns: {sorted(df.columns)}"
            )
        if min_area is not None:
            keep &= df["area"] >= min_area
        if max_area is not None:
            keep &= df["area"] <= max_area

    for prop, (lo, hi) in filters.items():
        if prop not in df.columns:
            raise ValueError(
                f"Unknown column '{prop}' in --filters. Available columns: "
                f"{sorted(df.columns)}"
            )
        if lo is not None:
            keep &= df[prop] >= lo
        if hi is not None:
            keep &= df[prop] <= hi

    return df[keep]


def build_relabel_maps(kept_ids_by_timepoint: dict[int | None, set]) -> dict[int | None, dict[int, int]]:
    """Map each kept old label id to a new, contiguous 1..N id, per timepoint."""
    return {
        t: {old_id: new_id for new_id, old_id in enumerate(sorted(ids), start=1)}
        for t, ids in kept_ids_by_timepoint.items()
    }


def remap_frame(frame: np.ndarray, mapping: dict[int, int]) -> np.ndarray:
    """Zero out every id not in `mapping`, and relabel the rest via `mapping`."""
    if not mapping:
        return np.zeros_like(frame)

    max_id = max(int(frame.max()) if frame.size else 0, max(mapping))
    lookup = np.zeros(max_id + 1, dtype=frame.dtype)
    for old_id, new_id in mapping.items():
        lookup[old_id] = new_id

    return lookup[frame]


def filter_label_level(
    level_data,
    t_index: int | None,
    relabel_maps: dict[int | None, dict[int, int]],
) -> da.Array:
    """Zero out every label id of `level_data` that isn't kept, and relabel
    the rest to a contiguous id space, per timepoint.

    Returned as a dask array explicitly chunked like `level_data`, so the
    written OME-Zarr level keeps the same chunksize as the input.
    """
    if t_index is None:
        array = np.asarray(level_data)
        filtered = remap_frame(array, relabel_maps[None]).astype(level_data.dtype)
    else:
        filtered = np.zeros(level_data.shape, dtype=level_data.dtype)
        for t in range(level_data.shape[t_index]):
            indexer = [slice(None)] * level_data.ndim
            indexer[t_index] = t
            frame = np.asarray(level_data[tuple(indexer)])
            filtered[tuple(indexer)] = remap_frame(frame, relabel_maps.get(t, {}))

    return da.from_array(filtered, chunks=level_data.chunksize)


def relabel_kept_rows(
    kept: pd.DataFrame,
    relabel_maps: dict[int | None, dict[int, int]],
    t_index: int | None,
) -> pd.DataFrame:
    kept = kept.copy()
    label_pos = list(kept.columns).index("label")

    if t_index is not None:
        new_labels = [
            relabel_maps[int(t)][old_id]
            for t, old_id in zip(kept["timepoint"], kept["label"])
        ]
    else:
        mapping = relabel_maps[None]
        new_labels = [mapping[old_id] for old_id in kept["label"]]

    kept["original_label"] = kept["label"]
    kept["label"] = new_labels

    # Put "original_label" right after "label", wherever "label" originally was.
    columns = [c for c in kept.columns if c != "original_label" and c != "label"]
    columns[label_pos:label_pos] = ["label", "original_label"]
    return kept[columns]


def main():
    args = parse_args()

    if args.output.exists():
        raise FileExistsError(f"Output already exists: {args.output}")

    filters = parse_filters(args.filters)

    df = pd.read_csv(args.csv)
    if "label" not in df.columns:
        raise ValueError(
            f"--csv {args.csv} has no 'label' column. Available columns: "
            f"{sorted(df.columns)}"
        )

    kept = filter_rows(df, args.min_area, args.max_area, filters)

    source = mm.ZarrManager(str(args.labels))
    metadata = copy.deepcopy(source.metadata)
    axes = metadata["axes"].lower()
    t_index = axes.index("t") if "t" in axes else None

    if t_index is not None:
        if "timepoint" not in kept.columns:
            raise ValueError(
                f"--labels {args.labels} has a time axis (axes={axes}) but "
                "--csv has no 'timepoint' column to group kept labels by."
            )
        kept_ids_by_timepoint = {
            int(t): set(group["label"].tolist())
            for t, group in kept.groupby("timepoint")
        }
    else:
        kept_ids_by_timepoint = {None: set(kept["label"].tolist())}

    relabel_maps = build_relabel_maps(kept_ids_by_timepoint)

    filtered_levels = [
        filter_label_level(level_data, t_index, relabel_maps)
        for level_data in source.data
    ]

    output = mm.ArrayManager(filtered_levels, metadata)
    output.to_zarr(str(args.output), overwrite=True)

    source.close()

    kept = relabel_kept_rows(kept, relabel_maps, t_index)
    kept.to_csv(args.output_csv, index=False)

    print(f"Kept {len(kept)}/{len(df)} label(s) across {len(kept_ids_by_timepoint)} timepoint(s)")
    print("Relabeled kept objects to a contiguous 1..N id space per timepoint")
    print(f"Wrote filtered labels to {args.output}")
    print(f"Wrote filtered measurements to {args.output_csv}")


if __name__ == "__main__":
    main()
