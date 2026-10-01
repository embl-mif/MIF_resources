#!/usr/bin/env python3
"""Preallocate an empty OME-Zarr skeleton.

Two modes:

1. --input-zarr given: inherit shape/axes/scale/dtype-adjacent metadata from
   an existing OME-Zarr store (any combination of NGFF 0.4/0.5, Zarr format
   2/3). Optionally override the channel count/labels (--num-channels/
   --channel-labels), restrict to a crop window per axis (--<axis>-min/
   --<axis>-max, in voxel coordinates at the finest *included* pyramid
   level), restrict the pyramid depth (--num-levels, --start-level), drop an
   axis from the output entirely (--drop-axes, e.g. 'c' for a label mask),
   and/or override the inherited chunk shape (--chunks, given in terms of
   the OUTPUT's own axes after any --drop-axes -- useful to align chunk
   boundaries with a tiling stride for safe concurrent writes into a shared
   merged store).

2. No --input-zarr: build a fresh, single-level skeleton directly from
   explicit --shape/--axes/--scale/--chunks/--translation geometry, with no
   source to inherit metadata from.
"""

import argparse
import math

import zarr
from ngio.common._pyramid import ImagePyramidBuilder

from ome_zarr_store import OMEZarrMetadata

CROP_AXES = ("t", "z", "y", "x")
SPATIAL_AXES = ("z", "y", "x")


def _axis_type(name: str) -> str:
    if name == "c":
        return "channel"
    if name == "t":
        return "time"
    return "space"


def _omero_channels(num_channels: int, channel_labels: list[str] | None, dtype: str) -> list[dict]:
    max_val = 255.0 if dtype == "uint8" else 1.0
    palette = ["FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "FFFFFF"]
    labels = channel_labels or [f"channel_{i}" for i in range(num_channels)]
    return [
        {
            "active": True,
            "coefficient": 1.0,
            "color": palette[i % len(palette)],
            "family": "linear",
            "inverted": False,
            "label": label_name,
            "window": {"end": max_val, "max": max_val, "min": 0.0, "start": 0.0},
        }
        for i, label_name in enumerate(labels)
    ]


def _rescale_crop_bound(value: float, base_scale: float, level_scale: float) -> float:
    """Rescale one crop endpoint given in the finest included level's voxel
    coordinates into another level's voxel coordinates, using the ratio of
    that level's own voxel scale to the base level's."""
    if level_scale == base_scale:
        return value
    return value * base_scale / level_scale


def build_skeleton(
    output_path: str,
    dtype: str,
    input_path: str | None = None,
    num_channels: int | None = None,
    channel_labels: list[str] | None = None,
    zarr_format: int | None = None,
    data_type: str | None = None,
    num_levels: int | None = None,
    crop: dict[str, tuple[int | None, int | None]] | None = None,
    start_level: int = 0,
    drop_axes: set[str] | None = None,
    chunks_override: list[int] | None = None,
    shape: list[int] | None = None,
    axis_names: list[str] | None = None,
    scale: list[float] | None = None,
    chunks: list[int] | None = None,
    translation: list[float] | None = None,
) -> None:
    crop = {axis: bounds for axis, bounds in (crop or {}).items() if bounds != (None, None)}
    drop_axes = drop_axes or set()

    if input_path is not None:
        _build_from_input(
            output_path, dtype, input_path, num_channels, channel_labels,
            zarr_format, data_type, num_levels, crop, start_level, drop_axes,
            chunks_override,
        )
        return

    if start_level:
        raise ValueError("--start-level requires --input-zarr.")
    if drop_axes:
        raise ValueError("--drop-axes requires --input-zarr.")
    if num_channels is not None:
        raise ValueError(
            "--num-channels requires --input-zarr -- with no source there's no default "
            "channel count to override; put the channel count in --shape directly."
        )
    if crop:
        raise ValueError(
            "--<axis>-min/--<axis>-max crop bounds require --input-zarr -- with no source "
            "there's no extent to crop from; state the cropped shape directly in --shape."
        )
    if num_levels is not None and num_levels != 1:
        raise ValueError(
            f"No --input-zarr given: output is always single-level, got --num-levels {num_levels}."
        )
    missing = [name for name, val in [("--shape", shape), ("--axes", axis_names),
                                       ("--scale", scale), ("--chunks", chunks)] if val is None]
    if missing:
        raise ValueError(f"No --input-zarr given: {', '.join(missing)} required.")

    _build_explicit(
        output_path, dtype, shape, axis_names, scale, chunks,
        translation, zarr_format, data_type, channel_labels,
    )


def _build_from_input(output_path: str, dtype: str, input_path: str,
                       num_channels: int | None, channel_labels: list[str] | None,
                       zarr_format: int | None, data_type: str | None,
                       num_levels: int | None,
                       crop: dict[str, tuple[int | None, int | None]],
                       start_level: int = 0,
                       drop_axes: set[str] | None = None,
                       chunks_override: list[int] | None = None) -> None:
    drop_axes = drop_axes or set()
    src_group = zarr.open(input_path, mode="r")
    md = OMEZarrMetadata(src_group)

    zattrs = dict(src_group.attrs)
    if data_type is not None:
        zattrs["data_type"] = data_type

    axis_names = tuple(axis["name"] for axis in md.axes)
    channel_axis_index = axis_names.index("c") if "c" in axis_names else None

    if num_channels is not None and channel_axis_index is None:
        raise ValueError(
            f"Input axes {axis_names} contain no 'c' axis; cannot place {num_channels} "
            "output channels."
        )
    if num_channels is not None and "c" in drop_axes:
        raise ValueError("--num-channels and --drop-axes c are mutually exclusive.")
    if channel_labels is not None and "c" in drop_axes:
        raise ValueError("--channel-labels and --drop-axes c are mutually exclusive.")
    for axis in crop:
        if axis not in axis_names:
            raise ValueError(
                f"Input axes {axis_names} contain no '{axis}' axis; cannot apply "
                f"--{axis}-min/--{axis}-max."
            )

    if zarr_format is None:
        zarr_format = md.zarr_format

    all_datasets = md.datasets  # [{"path", "scale", "translation"}, ...] -- already resolved
    if not (0 <= start_level < len(all_datasets)):
        raise ValueError(
            f"--start-level {start_level} out of range for the input's own "
            f"{len(all_datasets)} pyramid level(s)."
        )
    datasets = all_datasets[start_level:]

    if num_levels is not None:
        if num_levels < 1:
            raise ValueError(f"--num-levels must be >= 1, got {num_levels}")
        if num_levels > len(datasets):
            raise ValueError(
                f"--num-levels {num_levels} exceeds the {len(datasets)} pyramid level(s) "
                f"available from --start-level {start_level}."
            )
        datasets = datasets[:num_levels]
    read_level_paths = [ds["path"] for ds in datasets]
    write_level_paths = [str(i) for i in range(len(datasets))]

    base_array = src_group[read_level_paths[0]]
    base_scale = tuple(datasets[0]["scale"])

    base_crop: dict[str, tuple[int, int]] = {}
    for axis, (crop_min, crop_max) in crop.items():
        axis_index = axis_names.index(axis)
        extent = base_array.shape[axis_index]
        start = crop_min if crop_min is not None else 0
        end = crop_max if crop_max is not None else extent
        if not (0 <= start < end <= extent):
            raise ValueError(f"Bad {axis} range [{start}, {end}) for extent {extent}.")
        base_crop[axis] = (start, end)
        kind = "time" if axis == "t" else "spatial"
        print(f"[Prealloc] Restricting {axis} axis ({kind}) to [{start}, {end}) of {extent}")

    shapes = []
    scales = []
    translations = []
    for ds in datasets:
        src_array = src_group[ds["path"]]
        level_shape = list(src_array.shape)
        level_scale = ds["scale"]
        level_translation = list(ds["translation"])

        if num_channels is not None:
            level_shape[channel_axis_index] = num_channels

        for axis, (base_start, base_end) in base_crop.items():
            axis_index = axis_names.index(axis)
            start = _rescale_crop_bound(base_start, base_scale[axis_index], level_scale[axis_index])
            end = _rescale_crop_bound(base_end, base_scale[axis_index], level_scale[axis_index])
            start, end = int(math.floor(start)), min(int(math.ceil(end)), src_array.shape[axis_index])
            level_shape[axis_index] = end - start
            if axis in SPATIAL_AXES:
                level_translation[axis_index] += start * level_scale[axis_index]

        shapes.append(tuple(level_shape))
        scales.append(tuple(level_scale))
        translations.append(tuple(level_translation))

    if "t" in base_crop:
        zattrs["t_range"] = list(base_crop["t"])

    base_chunks = list(base_array.chunks)
    if num_channels is not None:
        base_chunks[channel_axis_index] = num_channels
    for axis in base_crop:
        axis_index = axis_names.index(axis)
        base_chunks[axis_index] = min(base_chunks[axis_index], shapes[0][axis_index])

    keep_indices = [i for i, name in enumerate(axis_names) if name not in drop_axes]
    output_axis_names = tuple(axis_names[i] for i in keep_indices)
    shapes = [tuple(shape[i] for i in keep_indices) for shape in shapes]
    scales = [tuple(scale[i] for i in keep_indices) for scale in scales]
    translations = [tuple(t[i] for i in keep_indices) for t in translations]
    base_chunks = [base_chunks[i] for i in keep_indices]
    if chunks_override is not None:
        if len(chunks_override) != len(output_axis_names):
            raise ValueError(
                f"--chunks has {len(chunks_override)} entries, but the output has "
                f"{len(output_axis_names)} axes {output_axis_names} (after any --drop-axes)."
            )
        base_chunks = list(chunks_override)

    builder = ImagePyramidBuilder.from_shapes(
        shapes=shapes,
        base_scale=scales,
        axes=output_axis_names,
        level_paths=write_level_paths,
        chunks=tuple(base_chunks),
        data_type=dtype,
        zarr_format=zarr_format,
    )

    dest_root = zarr.open(output_path, mode="w", zarr_format=zarr_format)
    builder.to_zarr(dest_root)

    ome_version = "0.5" if zarr_format == 3 else "0.4"

    new_datasets = []
    for path_name, level_scale, level_translation in zip(write_level_paths, scales, translations):
        transforms = [{"type": "scale", "scale": list(level_scale)}]
        if any(level_translation):
            transforms.append({"type": "translation", "translation": list(level_translation)})
        new_datasets.append({"path": path_name, "coordinateTransformations": transforms})

    # Carry each axis' unit over from the source so the output keeps the same physical units.
    source_units = {axis["name"]: axis.get("unit") for axis in md.axes}
    output_axes = []
    for name in output_axis_names:
        axis_entry = {"name": name, "type": _axis_type(name)}
        if source_units.get(name):
            axis_entry["unit"] = source_units[name]
        output_axes.append(axis_entry)

    new_ms = {
        "version": ome_version,
        "name": "/",
        "axes": output_axes,
        "datasets": new_datasets,
        "type": md.downsample_method,
    }
    new_ome_attrs = {"multiscales": [new_ms]}

    if "c" in drop_axes:
        pass  # no omero block for a channel-less output
    elif "c" in output_axis_names and (num_channels is not None or channel_labels is not None):
        channel_count = num_channels if num_channels is not None else shapes[0][output_axis_names.index("c")]
        new_ome_attrs["omero"] = {"channels": _omero_channels(channel_count, channel_labels, dtype)}

    # Always write multiscales/omero flat at the attrs root, regardless of
    # output zarr_format. This is NOT the real NGFF spec (which nests under
    # "ome" for zarr v3) -- it matches squirrel's own OMEZarrMetadata reader/
    # writer instead, which downstream tools in this pipeline (e.g.
    # populate_pyramid.py, via squirrel.OMEZarrStore) actually rely on:
    # squirrel's detect_ome_version() only ever checks attrs["multiscales"]
    # directly and never looks for a nested "ome" key, even for zarr v3.
    zattrs.pop("multiscales", None)
    zattrs.pop("omero", None)
    zattrs.pop("ome", None)
    zattrs.update(new_ome_attrs)

    dest_root.attrs.update(zattrs)

    for path_name, shape in zip(write_level_paths, shapes):
        print(f"[Prealloc] Created scale '{path_name}': shape {shape}")
    print(f"[Prealloc] Output Zarr format: v{zarr_format}")


def _build_explicit(output_path: str, dtype: str,
                     shape: list[int], axis_names: list[str],
                     scale: list[float], chunks: list[int],
                     translation: list[float] | None,
                     zarr_format: int, data_type: str | None,
                     channel_labels: list[str] | None) -> None:
    """No source zarr to inherit metadata from -- build a fresh, minimal
    skeleton straight from caller-specified geometry."""
    ndim = len(shape)
    if len(axis_names) != ndim:
        raise ValueError(f"--axes has {len(axis_names)} entries, --shape has {ndim}")
    if len(scale) != ndim:
        raise ValueError(f"--scale has {len(scale)} entries, --shape has {ndim}")
    if len(chunks) != ndim:
        raise ValueError(f"--chunks has {len(chunks)} entries, --shape has {ndim}")

    translation = translation if translation is not None else [0.0] * ndim
    if len(translation) != ndim:
        raise ValueError(f"--translation has {len(translation)} entries, --shape has {ndim}")

    builder = ImagePyramidBuilder.from_shapes(
        shapes=[tuple(shape)],
        base_scale=[tuple(scale)],
        axes=tuple(axis_names),
        level_paths=["0"],
        chunks=tuple(chunks),
        data_type=dtype,
        zarr_format=zarr_format,
    )

    dest_root = zarr.open(output_path, mode="w", zarr_format=zarr_format)
    builder.to_zarr(dest_root)

    ome_version = "0.5" if zarr_format == 3 else "0.4"

    zattrs = {
        "multiscales": [{
            "version": ome_version,
            "name": "/",
            "axes": [{"name": name, "type": _axis_type(name)} for name in axis_names],
            "datasets": [{
                "path": "0",
                "coordinateTransformations": [
                    {"type": "scale", "scale": list(scale)},
                    {"type": "translation", "translation": list(translation)},
                ],
            }],
            "type": "identity",
        }]
    }
    if data_type is not None:
        zattrs["data_type"] = data_type

    if "c" in axis_names and channel_labels is not None:
        c_axis = axis_names.index("c")
        zattrs["omero"] = {"channels": _omero_channels(shape[c_axis], channel_labels, dtype)}

    dest_root.attrs.update(zattrs)

    print(f"[Prealloc] Created single-level skeleton: shape {tuple(shape)}")
    print(f"[Prealloc] Output Zarr format: v{zarr_format}")


def _parse_csv_ints(text: str) -> list[int]:
    return [int(v) for v in text.split(",")]


def _parse_csv_floats(text: str) -> list[float]:
    return [float(v) for v in text.split(",")]


def _parse_csv_strs(text: str) -> list[str]:
    return [v for v in text.split(",")]


def main():
    parser = argparse.ArgumentParser(
        description="Preallocate an empty OME-Zarr skeleton -- either derived from an "
                     "existing input zarr's own metadata, or from explicit geometry."
    )
    parser.add_argument("--input-zarr", default=None,
                         help="Derive geometry and inherit metadata from this existing "
                              "zarr. Mutually exclusive with --shape/--axes/--scale/"
                              "--chunks, which are required if this is omitted.")
    parser.add_argument("--output-zarr", required=True)
    parser.add_argument("--dtype", choices=["float32", "uint8", "int32", "uint16", "uint32"],
                         default="float32")
    parser.add_argument("--num-channels", type=int, default=None,
                         help="Requires --input-zarr. Number of output channels, "
                              "overriding the source's own channel count.")
    parser.add_argument("--channel-labels", type=str, default=None,
                         help="Comma-separated channel names. With --input-zarr, matches "
                              "--num-channels in length if given, else the source's own "
                              "channel count. Without --input-zarr, matches --shape's "
                              "channel-axis size.")
    parser.add_argument("--zarr-format", type=int, choices=[2, 3], default=None,
                         help="Output Zarr storage format. With --input-zarr, defaults to "
                              "matching its format. REQUIRED without --input-zarr.")
    parser.add_argument("--data-type", type=str, default=None,
                         choices=["intensity", "label", "probability"])
    parser.add_argument("--num-levels", type=int, default=None,
                         help="Requires --input-zarr. Restrict to the first N pyramid "
                              "levels (counted from --start-level, if given). Without "
                              "--input-zarr, output is always single-level.")
    parser.add_argument("--start-level", type=int, default=0,
                         help="Requires --input-zarr. Treat this index into the source's "
                              "own pyramid as the output's level 0, dropping everything "
                              "finer. Destination levels are always renumbered "
                              "'0','1',... regardless of this value; reads still use the "
                              "source's own path names.")
    parser.add_argument("--drop-axes", type=str, default=None,
                         help="Requires --input-zarr. Comma-separated axis names to "
                              "remove from the output geometry entirely (e.g. 'c' for a "
                              "label mask that never has a channel axis). Lenient about "
                              "axes the source doesn't have.")
    for axis in CROP_AXES:
        parser.add_argument(f"--{axis}-min", type=int, default=None,
                             help=f"Requires --input-zarr. Crop the {axis} axis to start here "
                                  f"(voxel coords at the finest included level).")
        parser.add_argument(f"--{axis}-max", type=int, default=None,
                             help=f"Requires --input-zarr. Crop the {axis} axis to end here "
                                  f"(exclusive, voxel coords at the finest included level).")
    parser.add_argument("--shape", type=str, default=None,
                         help="No --input-zarr: comma-separated ints, e.g. '10,1024,1024'.")
    parser.add_argument("--axes", type=str, default=None,
                         help="No --input-zarr: comma-separated axis names matching "
                              "--shape, e.g. 'z,y,x'.")
    parser.add_argument("--scale", type=str, default=None,
                         help="No --input-zarr: comma-separated per-axis voxel size, "
                              "matching --shape.")
    parser.add_argument("--translation", type=str, default=None,
                         help="No --input-zarr: comma-separated per-axis physical-space "
                              "origin, matching --shape. Defaults to all zeros if omitted.")
    parser.add_argument("--chunks", type=str, default=None,
                         help="Comma-separated per-axis chunk size. With --input-zarr, "
                              "optional -- overrides the source's own inherited chunk "
                              "shape, given in terms of the OUTPUT's axes (after any "
                              "--drop-axes), e.g. to align with a tiling stride. Without "
                              "--input-zarr, required, matching --shape.")
    args = parser.parse_args()

    channel_labels = _parse_csv_strs(args.channel_labels) if args.channel_labels else None
    drop_axes = set(_parse_csv_strs(args.drop_axes)) if args.drop_axes else set()
    crop = {axis: (getattr(args, f"{axis}_min"), getattr(args, f"{axis}_max")) for axis in CROP_AXES}
    crop_given = {axis: bounds for axis, bounds in crop.items() if bounds != (None, None)}
    explicit_mode_args_given = any([args.shape, args.axes, args.scale])

    if args.input_zarr and explicit_mode_args_given:
        parser.error("--input-zarr and --shape/--axes/--scale/--chunks are mutually "
                     "exclusive.")
    if not args.input_zarr and not explicit_mode_args_given:
        parser.error("Provide either --input-zarr or --shape/--axes/--scale/--chunks.")

    if args.input_zarr:
        if channel_labels and args.num_channels is not None and len(channel_labels) != args.num_channels:
            parser.error(f"--channel-labels has {len(channel_labels)} entries, "
                         f"but --num-channels is {args.num_channels}")
        if args.num_levels is not None and args.num_levels < 1:
            parser.error(f"--num-levels must be >= 1, got {args.num_levels}")
        if args.start_level < 0:
            parser.error(f"--start-level must be >= 0, got {args.start_level}")

        build_skeleton(
            output_path=args.output_zarr, dtype=args.dtype, input_path=args.input_zarr,
            num_channels=args.num_channels, channel_labels=channel_labels,
            zarr_format=args.zarr_format, data_type=args.data_type,
            num_levels=args.num_levels, crop=crop_given,
            start_level=args.start_level, drop_axes=drop_axes,
            chunks_override=_parse_csv_ints(args.chunks) if args.chunks else None,
        )
    else:
        if not (args.shape and args.axes and args.scale and args.chunks):
            parser.error("Requires --shape, --axes, --scale, and --chunks all set.")
        if args.zarr_format is None:
            parser.error("--zarr-format is required without --input-zarr.")
        if args.num_levels is not None and args.num_levels != 1:
            parser.error(f"Output is always single-level without --input-zarr; "
                         f"--num-levels must be 1 if given, got {args.num_levels}.")
        if args.num_channels is not None:
            parser.error("--num-channels requires --input-zarr -- channel count comes "
                         "from --shape directly otherwise.")
        if crop_given:
            parser.error("--<axis>-min/--<axis>-max require --input-zarr -- state the "
                         "cropped shape directly in --shape otherwise.")
        if args.start_level:
            parser.error("--start-level requires --input-zarr.")
        if drop_axes:
            parser.error("--drop-axes requires --input-zarr.")

        build_skeleton(
            output_path=args.output_zarr, dtype=args.dtype,
            shape=_parse_csv_ints(args.shape), axis_names=_parse_csv_strs(args.axes),
            scale=_parse_csv_floats(args.scale), chunks=_parse_csv_ints(args.chunks),
            translation=_parse_csv_floats(args.translation) if args.translation else None,
            zarr_format=args.zarr_format, data_type=args.data_type,
            channel_labels=channel_labels,
        )

    print(f"[Prealloc] Done: {args.output_zarr}")


if __name__ == "__main__":
    main()
