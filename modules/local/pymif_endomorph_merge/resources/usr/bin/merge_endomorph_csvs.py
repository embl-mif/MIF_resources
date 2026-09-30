#!/usr/bin/env python3

"""
Concatenate the per-dataset endo_morph CSVs (per-cell and per-image summary)
into one CSV each. Rows keep their `dataset_id` column; columns are the union
of the inputs (missing values left empty).
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", nargs="+", required=True, type=Path)
    parser.add_argument("--summaries", nargs="+", required=True, type=Path)
    parser.add_argument("--output-cells", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def merge(paths: list[Path], output: Path) -> None:
    # Sorted so the merged rows do not depend on task completion order.
    frames = [pd.read_csv(p) for p in sorted(paths)]
    frames = [f for f in frames if len(f)] or frames[:1]
    merged = pd.concat(frames, ignore_index=True, sort=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    print(f"Merged {len(paths)} file(s), {len(merged)} row(s) -> {output}")


def main():
    args = parse_args()
    merge(args.cells, args.output_cells)
    merge(args.summaries, args.output_summary)


if __name__ == "__main__":
    main()
