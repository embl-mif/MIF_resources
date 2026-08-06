#!/usr/bin/env python3

"""Construct regular tiles for an OME-Zarr image and write tile bounds to CSV.

The output bounds are always written in voxel units as half-open intervals
[min, max) for each spatial axis.

Rows are emitted in the canonical Cellpose table shape to avoid downstream
schema adaptation. The dataset_id is unique per tile and formed as
<source_dataset_id>__<tile_id>.
"""

import argparse
import csv
import math
from typing import Any

import ngff_zarr as nz


def _normalize_bounds_arg(
    min_value: float | None,
    max_value: float | None,
    axis_name: str,
) -> tuple[float | None, float | None] | None:
    if min_value is None and max_value is None:
        return None

    start = float(min_value) if min_value is not None else None
    end = float(max_value) if max_value is not None else None
    if start is not None and end is not None and end <= start:
        raise ValueError(f"{axis_name}-bounds are invalid: min={start} must be < max={end}")
    return start, end


def _resolve_spatial_roi(
    data: Any,
    dims: list[str],
    scales: dict[str, float],
    translation: dict[str, float],
    x_min: float | None,
    x_max: float | None,
    y_min: float | None,
    y_max: float | None,
    z_min: float | None,
    z_max: float | None,
    use_physical_units: bool,
) -> dict[str, tuple[int, int]]:
    x_bounds = _normalize_bounds_arg(x_min, x_max, "x")
    y_bounds = _normalize_bounds_arg(y_min, y_max, "y")
    z_bounds = _normalize_bounds_arg(z_min, z_max, "z")

    if "x" not in dims or "y" not in dims:
        raise ValueError(f"Input dims must contain x and y axes, got {dims}")

    has_z_axis = "z" in dims
    if z_bounds is not None and not has_z_axis:
        raise ValueError("Input has no z dimension, but --z-min/--z-max was provided")

    resolved_bounds: dict[str, tuple[int, int]] = {}
    for dim in ("x", "y", "z"):
        if dim not in dims:
            continue

        axis = dims.index(dim)
        axis_size = int(data.shape[axis])
        scale = float(scales.get(dim, 1.0))
        if scale <= 0.0:
            raise ValueError(f"{dim} scale must be > 0, got {scale}")
        offset = float(translation.get(dim, 0.0))

        if dim == "x":
            interval = x_bounds
        elif dim == "y":
            interval = y_bounds
        else:
            interval = z_bounds

        if interval is None:
            start, end = 0, axis_size
        else:
            start_input, end_input = interval
            if use_physical_units:
                if start_input is None:
                    start = 0
                else:
                    relative_start = (start_input - offset) / scale
                    start = int(math.floor(relative_start + 1e-9))

                if end_input is None:
                    end = axis_size
                else:
                    relative_end = (end_input - offset) / scale
                    end = int(math.ceil(relative_end - 1e-9))

                if (start_input is not None and start < 0) or (end_input is not None and end > axis_size):
                    raise ValueError(
                        f"{dim}-bounds [{start_input}, {end_input}) are outside the image extent "
                        f"[{offset}, {offset + (axis_size * scale)}) at the selected resolution level"
                    )
            else:
                if start_input is None:
                    start = 0
                else:
                    start = int(math.floor(start_input + 1e-9))

                if end_input is None:
                    end = axis_size
                else:
                    end = int(math.ceil(end_input - 1e-9))

                if (start_input is not None and start < 0) or (end_input is not None and end > axis_size):
                    raise ValueError(
                        f"{dim}-bounds [{start_input}, {end_input}) in voxel units are outside the image extent "
                        f"[0, {axis_size}) at the selected resolution level"
                    )

            if end <= start:
                raise ValueError(
                    f"{dim}-bounds produce an empty ROI after conversion: start={start}, end={end}"
                )

        resolved_bounds[dim] = (start, end)

    return resolved_bounds


def _to_voxel_size(size: float, scale: float, axis_name: str, use_physical_units: bool) -> int:
    if size <= 0.0:
        raise ValueError(f"tile size for axis '{axis_name}' must be > 0, got {size}")
    if use_physical_units:
        voxels = int(math.ceil((size / scale) - 1e-9))
    else:
        voxels = int(math.ceil(size - 1e-9))
    return max(1, voxels)


def _to_voxel_overlap(overlap: float, scale: float, axis_name: str, use_physical_units: bool) -> int:
    if overlap < 0.0:
        raise ValueError(f"tile overlap for axis '{axis_name}' must be >= 0, got {overlap}")
    if use_physical_units:
        voxels = int(math.ceil((overlap / scale) - 1e-9))
    else:
        voxels = int(math.ceil(overlap - 1e-9))
    return max(0, voxels)


def _axis_tiles(start: int, end: int, tile_size: int, overlap: int) -> list[tuple[int, int]]:
    if tile_size <= overlap:
        raise ValueError(
            f"tile_size ({tile_size}) must be greater than overlap ({overlap}) for each tiled axis"
        )

    axis_length = end - start
    if axis_length <= 0:
        raise ValueError(f"invalid axis bounds [{start}, {end})")

    step = tile_size - overlap
    tiles: list[tuple[int, int]] = []
    current_start = start

    while current_start < end:
        # Keep boundary tiles at full size even when they extend beyond the ROI/image.
        current_end = current_start + tile_size
        tiles.append((current_start, current_end))
        if current_start + step >= end:
            break
        current_start += step

    return tiles


def construct_tiles(
    source_dataset_id: str,
    input_image_path: str,
    output_csv: str,
    tile_size_x: float,
    tile_size_y: float,
    tile_size_z: float | None,
    tile_overlap: float,
    resolution_level: int,
    x_min: float | None,
    x_max: float | None,
    y_min: float | None,
    y_max: float | None,
    z_min: float | None,
    z_max: float | None,
    use_physical_units: bool,
) -> None:
    multiscales = nz.from_ngff_zarr(input_image_path)
    if resolution_level < 0:
        raise ValueError(f"--resolution-level must be >= 0, got {resolution_level}")
    if resolution_level >= len(multiscales.images):
        raise ValueError(
            f"--resolution-level {resolution_level} is out of range for input with {len(multiscales.images)} level(s)"
        )

    ngff_image = multiscales.images[resolution_level]
    data = ngff_image.data
    dims = list(ngff_image.dims)
    scales = {dim: float(ngff_image.scale.get(dim, 1.0)) for dim in dims}
    translation = {dim: float((getattr(ngff_image, "translation", {}) or {}).get(dim, 0.0)) for dim in dims}

    roi_bounds = _resolve_spatial_roi(
        data=data,
        dims=dims,
        scales=scales,
        translation=translation,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        z_min=z_min,
        z_max=z_max,
        use_physical_units=use_physical_units,
    )

    if "z" in dims and tile_size_z is None:
        raise ValueError("Input has z dimension; --tile-size-z must be provided")

    tile_sizes_voxel = {
        "x": _to_voxel_size(tile_size_x, scales["x"], "x", use_physical_units),
        "y": _to_voxel_size(tile_size_y, scales["y"], "y", use_physical_units),
    }
    if "z" in roi_bounds:
        if tile_size_z is None:
            raise ValueError("--tile-size-z is required when tiling a z axis")
        tile_sizes_voxel["z"] = _to_voxel_size(tile_size_z, scales["z"], "z", use_physical_units)

    overlap_voxel = {
        axis: _to_voxel_overlap(tile_overlap, scales[axis], axis, use_physical_units)
        for axis in roi_bounds
    }

    axis_tiles = {
        axis: _axis_tiles(
            start=roi_bounds[axis][0],
            end=roi_bounds[axis][1],
            tile_size=tile_sizes_voxel[axis],
            overlap=overlap_voxel[axis],
        )
        for axis in roi_bounds
    }

    rows: list[dict[str, Any]] = []
    tile_index = 0

    z_tiles = axis_tiles.get("z", [(0, 1)])
    for z_bounds_vox in z_tiles:
        for y_bounds_vox in axis_tiles["y"]:
            for x_bounds_vox in axis_tiles["x"]:
                tile_index += 1
                row: dict[str, Any] = {
                    "tile_id": f"tile_{tile_index:06d}",
                    "source_dataset_id": source_dataset_id,
                    "dataset_id": "",
                    "input_uri": input_image_path,
                    "x_min": x_bounds_vox[0],
                    "x_max": x_bounds_vox[1],
                    "y_min": y_bounds_vox[0],
                    "y_max": y_bounds_vox[1],
                    "z_min": "",
                    "z_max": "",
                    "resolution_level": int(resolution_level),
                }
                row["dataset_id"] = f"{source_dataset_id}__{row['tile_id']}"
                if "z" in roi_bounds:
                    row["z_min"] = z_bounds_vox[0]
                    row["z_max"] = z_bounds_vox[1]
                rows.append(row)

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_id",
                "source_dataset_id",
                "tile_id",
                "input_uri",
                "x_min",
                "x_max",
                "y_min",
                "y_max",
                "z_min",
                "z_max",
                "resolution_level",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} tile(s) to {output_csv}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()

    parser.add_argument("--source-dataset-id", type=str, required=True)
    parser.add_argument("--input-image-path", type=str, required=True)
    parser.add_argument("--output-csv", type=str, required=True)

    parser.add_argument("--tile-size-x", type=float, required=True)
    parser.add_argument("--tile-size-y", type=float, required=True)
    parser.add_argument("--tile-size-z", type=float, default=None)

    parser.add_argument("--tile-overlap", type=float, default=0.0)
    parser.add_argument("--resolution-level", type=int, default=0)

    parser.add_argument("--x-min", type=float, default=None)
    parser.add_argument("--x-max", type=float, default=None)
    parser.add_argument("--y-min", type=float, default=None)
    parser.add_argument("--y-max", type=float, default=None)
    parser.add_argument("--z-min", type=float, default=None)
    parser.add_argument("--z-max", type=float, default=None)

    parser.add_argument(
        "--use-physical-units",
        action="store_true",
        help="Interpret tile sizes, overlap, and ROI bounds in physical units.",
    )

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    construct_tiles(
        source_dataset_id=args.source_dataset_id,
        input_image_path=args.input_image_path,
        output_csv=args.output_csv,
        tile_size_x=args.tile_size_x,
        tile_size_y=args.tile_size_y,
        tile_size_z=args.tile_size_z,
        tile_overlap=args.tile_overlap,
        resolution_level=args.resolution_level,
        x_min=args.x_min,
        x_max=args.x_max,
        y_min=args.y_min,
        y_max=args.y_max,
        z_min=args.z_min,
        z_max=args.z_max,
        use_physical_units=args.use_physical_units,
    )


if __name__ == "__main__":
    main()

