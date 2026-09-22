#!/usr/bin/env python3

"""Compute the --start-level/--<axis>-min/--<axis>-max/--chunks flags that
let preallocate_zarr.py build the merged-segmentation OME-Zarr store as a
real crop of the original image's own pyramid.

Prints (to stdout, one line, space-separated) rather than building anything
itself, so it can be dropped into a Nextflow process's script block via
command substitution ahead of the actual preallocate_zarr.py call:

    params=$(calculate_zarr_params.py --tiles-csv "${tiles_csv}")
    preallocate_zarr.py --input-zarr "${source_zarr}" \
        --output-zarr "${output_zarr}" --drop-axes c \
        --dtype uint32 --data-type label $params
"""

import argparse

from zarr_tile_utils import _read_tile_rows, compute_output_geometry


def calculate_params(tiles_csv: str) -> dict:
    rows = _read_tile_rows(tiles_csv)
    if not rows:
        raise ValueError(f"{tiles_csv} has no tile rows")

    # resolution_level is a plain CSV column, constant across every row of one
    # source image (written once per image by tile_constructor.py) -- read it
    # here rather than threading it in separately as another CLI argument.
    resolution_level = int(rows[0]["resolution_level"])

    axis_order, shape, chunks, global_min = compute_output_geometry(rows)
    bounds = {
        axis: (global_min[i], global_min[i] + shape[i])
        for i, axis in enumerate(axis_order)
    }

    return {
        "resolution_level": resolution_level,
        "axis_order": axis_order,
        "bounds": bounds,
        "chunks": dict(zip(axis_order, chunks)),
    }


def format_as_flags(params: dict) -> str:
    """Render computed params as preallocate_zarr.py's own explicit-geometry
    flags -- comma-joined per flag, space-separated between flags, safe for
    a shell to word-split via `$(...)` command substitution."""
    flags = [f"--start-level {params['resolution_level']}"]
    for axis in params["axis_order"]:
        start, end = params["bounds"][axis]
        flags.append(f"--{axis}-min {start}")
        flags.append(f"--{axis}-max {end}")
        flags.append(f"--chunks {','.join(str(params['chunks'][axis]) for axis in params['axis_order'])}")
        return " ".join(flags)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute the merged-segmentation OME-Zarr store's --input-zarr-mode "
                     "flags (--start-level / crop bounds / --chunks) from its tile CSV's "
                     "own ROI bounds and resolution_level."
    )

    parser.add_argument("--tiles-csv", type=str, required=True)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    params = calculate_params(args.tiles_csv)
    print(format_as_flags(params))


if __name__ == "__main__":
    main()
