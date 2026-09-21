#!/usr/bin/env python3

"""Find forward spatial neighbors (right/bottom/behind) for every tile of one
source image and emit paths to their segmented label zarrs.

Pure metadata operation over `tiles_and_segmentations.csv` -- no zarr/image
I/O. Tile B is a forward neighbor of tile A along one axis if B starts inside
A's overlap band on that axis, and A and B match exactly on every other axis
(TILE_CONSTRUCTOR tiles each axis independently, so every tile in the same
row/column shares identical bounds on the other axes -- a mere overlap test
there would also match diagonal tiles, since adjacent rows always overlap by
`tile_overlap` near their shared border). Restricting to forward neighbors
(A's coordinate <= B's) means each adjacent pair is matched exactly once,
from the lower-coordinate tile's side.

tile_overlap_x/y/z are read per row from the input CSV rather than taken as a
CLI argument: TILE_CONSTRUCTOR already converts overlap to voxel units
per-axis (it can differ between axes under anisotropic physical-unit
scaling), so re-deriving or re-converting it here would risk disagreeing with
that conversion.
"""

import argparse
import csv
from typing import Any


def _to_int(value: str | None) -> int | None:
    text = value.strip() if value is not None else ""
    return int(text) if text else None


def _load_rows(tiles_csv: str) -> list[dict[str, Any]]:
    with open(tiles_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [
            {
                "source_dataset_id": row["source_dataset_id"],
                "dataset_id": row["dataset_id"],
                "tile_index": _to_int(row["tile_index"]),
                "segmented_tile_path": row["segmented_tile_path"],
                "x_min": _to_int(row["x_min"]),
                "x_max": _to_int(row["x_max"]),
                "y_min": _to_int(row["y_min"]),
                "y_max": _to_int(row["y_max"]),
                "z_min": _to_int(row["z_min"]),
                "z_max": _to_int(row["z_max"]),
                "tile_overlap_x": _to_int(row["tile_overlap_x"]),
                "tile_overlap_y": _to_int(row["tile_overlap_y"]),
                "tile_overlap_z": _to_int(row["tile_overlap_z"]),
            }
            for row in reader
        ]

    if not rows:
        raise ValueError(f"No rows found in {tiles_csv}")
    return rows


def _same_range(a_min: int, a_max: int, b_min: int, b_max: int) -> bool:
    return a_min == b_min and a_max == b_max


def _touches_forward(a_max: int, b_min: int, overlap: int) -> bool:
    return b_min >= a_max - overlap and b_min <= a_max


def find_neighbors(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    has_z = rows[0]["z_min"] is not None
    if any((row["z_min"] is not None) != has_z for row in rows):
        raise ValueError("Inconsistent z axis presence across tiles of the same source image")

    results: list[dict[str, Any]] = []
    for self_row in rows:
        right = bottom = behind = None

        for other in rows:
            if other is self_row:
                continue

            same_y = _same_range(self_row["y_min"], self_row["y_max"], other["y_min"], other["y_max"])
            same_x = _same_range(self_row["x_min"], self_row["x_max"], other["x_min"], other["x_max"])
            same_z = (not has_z) or _same_range(
                self_row["z_min"], self_row["z_max"], other["z_min"], other["z_max"]
            )

            if same_y and same_z and _touches_forward(self_row["x_max"], other["x_min"], self_row["tile_overlap_x"]):
                if right is not None:
                    raise ValueError(f"Multiple right neighbors found for tile {self_row['dataset_id']}")
                right = other

            if same_x and same_z and _touches_forward(self_row["y_max"], other["y_min"], self_row["tile_overlap_y"]):
                if bottom is not None:
                    raise ValueError(f"Multiple bottom neighbors found for tile {self_row['dataset_id']}")
                bottom = other

            if has_z and same_x and same_y and _touches_forward(
                self_row["z_max"], other["z_min"], self_row["tile_overlap_z"]
            ):
                if behind is not None:
                    raise ValueError(f"Multiple behind neighbors found for tile {self_row['dataset_id']}")
                behind = other

        results.append(
            {
                "source_dataset_id": self_row["source_dataset_id"],
                "dataset_id": self_row["dataset_id"],
                "segmented_tile_path": self_row["segmented_tile_path"],
                "right_dataset_id": right["dataset_id"] if right else "",
                "right_tile_index": right["tile_index"] if right else "",
                "right_segmented_tile_path": right["segmented_tile_path"] if right else "",
                "bottom_dataset_id": bottom["dataset_id"] if bottom else "",
                "bottom_tile_index": bottom["tile_index"] if bottom else "",
                "bottom_segmented_tile_path": bottom["segmented_tile_path"] if bottom else "",
                "behind_dataset_id": behind["dataset_id"] if behind else "",
                "behind_tile_index": behind["tile_index"] if behind else "",
                "behind_segmented_tile_path": behind["segmented_tile_path"] if behind else "",
            }
        )

    return results


def write_neighbors_csv(neighbor_rows: list[dict[str, Any]], output_csv: str) -> None:
    fieldnames = [
        "source_dataset_id",
        "dataset_id",
        "segmented_tile_path",
        "right_dataset_id",
        "right_tile_index",
        "right_segmented_tile_path",
        "bottom_dataset_id",
        "bottom_tile_index",
        "bottom_segmented_tile_path",
        "behind_dataset_id",
        "behind_tile_index",
        "behind_segmented_tile_path",
    ]
    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(neighbor_rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tiles-csv", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    rows = _load_rows(args.tiles_csv)
    neighbor_rows = find_neighbors(rows)
    write_neighbors_csv(neighbor_rows, args.output)
    print(f"Wrote {len(neighbor_rows)} tile neighbor row(s) to {args.output}")


if __name__ == "__main__":
    main()
