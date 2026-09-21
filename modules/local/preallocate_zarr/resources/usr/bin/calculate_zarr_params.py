#!/usr/bin/env python3

"""Compute the shape/axes/scale/chunks/translation/zarr-format for the
merged-segmentation OME-Zarr store of one source image, from its tile CSV's
own ROI bounds, and print them as ready-to-use preallocate_zarr.py flags.

Prints (to stdout, one line, space-separated) rather than building anything
itself, so it can be dropped into a Nextflow process's script block via
command substitution ahead of the actual preallocate_zarr.py call:

    params=$(calculate_zarr_params.py --tiles-csv "${tiles_csv}" \\
        --input-image-path "${input_image_path}")
    preallocate_zarr.py --output-zarr "${output_zarr}" --dtype uint32 \\
        --data-type label $params

That keeps `build_skeleton` in preallocate_zarr.py the only code path in
this module that ever actually creates a store -- this script only figures
out what to tell it. Lives alongside preallocate_zarr.py in the same
PREALLOCATE_ZARR module/bin; the old bespoke OMEZarrStore-based
preallocate_empty_zarr.py script/module is retired entirely.

Shape is the union (min/max) of every tile's ROI bounds in
`tiles_and_segmentations.csv` -- Module 1 already resolved this, it's never
re-derived from the original image. Per-axis chunk size is the tiling
stride (`tile_size - tile_overlap`), constant across every tile of one
source image (both are written once per source image by
`tile_constructor.py`, not re-derived per tile). Choosing the stride as the
chunk size is what makes every tile's later core-region write (Module 6b,
`write_relabeled_tile.py`) land exactly on the storage grid, with the
boundary/merged tile exempted by `check_grid_alignment`'s own "reaches the
true dataset edge" rule.

The output's OME-NGFF version and Zarr format match the *original* input
image's, detected via `zarr_tile_utils._read_input_calibration` (itself
built on `ome_zarr_store.detect_ome_version`/`detect_zarr_format`) rather
than hardcoded -- mirroring `bin/cellpose_sam_zarr.py`'s own
`detect_ome_zarr_version` precedent of propagating the input's version to
its output. `--zarr-format` overrides this detection, same as
preallocate_zarr.py's own --input-zarr path.
"""

import argparse

from zarr_tile_utils import _read_input_calibration, _read_tile_rows, compute_output_geometry


def calculate_params(tiles_csv: str, input_image_path: str,
                      zarr_format_override: int | None) -> dict:
    rows = _read_tile_rows(tiles_csv)
    if not rows:
        raise ValueError(f"{tiles_csv} has no tile rows")

    # resolution_level is a plain CSV column, constant across every row of one
    # source image (written once per image by tile_constructor.py) -- read it
    # here rather than threading it in separately as another CLI argument.
    resolution_level = int(rows[0]["resolution_level"])

    axis_order, shape, chunks, global_min = compute_output_geometry(rows)
    ome_version, detected_zarr_format, resolution = _read_input_calibration(
        input_image_path, resolution_level, axis_order
    )
    zarr_format = zarr_format_override if zarr_format_override is not None else detected_zarr_format

    # NGFF's "translation" coordinateTransformation is applied after "scale",
    # so it's expressed in the same physical units -- the ROI's voxel origin
    # (`global_min`) times the voxel size (`resolution`), not the raw voxel
    # offset itself. Same convention preallocate_zarr.py's own --input-zarr
    # crop handling uses for spatial axes.
    translation = tuple(origin * scale for origin, scale in zip(global_min, resolution))

    return {
        "shape": shape,
        "axes": axis_order,
        "scale": resolution,
        "chunks": chunks,
        "translation": translation,
        "zarr_format": zarr_format,
        "ome_version": ome_version,
    }


def format_as_flags(params: dict) -> str:
    """Render computed params as preallocate_zarr.py's own explicit-geometry
    flags -- comma-joined per flag, space-separated between flags, safe for
    a shell to word-split via `$(...)` command substitution."""
    return " ".join([
        f"--shape {','.join(str(v) for v in params['shape'])}",
        f"--axes {','.join(params['axes'])}",
        f"--scale {','.join(str(v) for v in params['scale'])}",
        f"--chunks {','.join(str(v) for v in params['chunks'])}",
        f"--translation {','.join(str(v) for v in params['translation'])}",
        f"--zarr-format {params['zarr_format']}",
    ])


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute the merged-segmentation OME-Zarr store's geometry from its "
                     "tile CSV's own ROI bounds, and print it as preallocate_zarr.py flags."
    )

    parser.add_argument("--tiles-csv", type=str, required=True)
    parser.add_argument("--input-image-path", type=str, required=True)
    parser.add_argument("--zarr-format", type=int, choices=[2, 3], default=None,
                         help="Override the detected output Zarr storage format. Defaults "
                              "to matching the input image's own format (detected, not "
                              "assumed).")

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    params = calculate_params(args.tiles_csv, args.input_image_path, args.zarr_format)
    print(format_as_flags(params))


if __name__ == "__main__":
    main()
