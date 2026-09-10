#!/usr/bin/env python3
"""Preallocate an empty OME-Zarr skeleton, in one of two modes:

1. Input-derived (--input-zarr): copy an existing source zarr's own
   metadata (axes, scale, omero channels, data_type), optionally remapping
   the channel axis or restricting the t axis / pyramid level count. Used
   for Ilastik probability-map stores and Ultrack inputs -- both are
   derivatives of one specific source image, so inheriting its metadata
   (minus targeted overrides) is correct.

2. Explicit-geometry (--shape/--axes/--scale/--chunks): build a completely
   fresh, minimal skeleton from caller-specified geometry, with no source
   zarr to copy metadata from. Always single-resolution (no pyramid --
   without a known input to derive a downsample-factor progression from,
   there's no well-defined way to build further levels). This is the mode
   a shared merge target that many independent tasks each write one
   disjoint region into needs -- e.g. PREALLOCATE_EMPTY_ZARR's stitched
   multi-tile label store, which cannot use mode 1 because there's no
   single source image whose metadata should be copied wholesale onto a
   uint32 label array.

"""

import argparse

import zarr
from ngio.common._pyramid import ImagePyramidBuilder
from ome_zarr_store import detect_ome_version, detect_zarr_format


def build_skeleton_from_input(input_path: str, output_path: str, dtype: str,
                               num_channels: int | None, channel_labels: list[str] | None,
                               zarr_format: int | None, data_type: str | None,
                               t_min: int | None, t_max: int | None,
                               num_levels: int | None = None):
    """Mode 1: derive geometry AND inherit metadata from an existing input zarr."""
    src_group = zarr.open(input_path, mode="r")
    zattrs = dict(src_group.attrs)
    multiscales = zattrs.get("multiscales", [])
    if not multiscales:
        raise ValueError(f"Input Zarr '{input_path}' lacks 'multiscales' metadata.")
    if data_type is not None:
        zattrs["data_type"] = data_type

    axes_meta = multiscales[0].get("axes", [])
    axis_names = tuple(
        (axis.get("name") if isinstance(axis, dict) else axis) for axis in axes_meta
    )
    channel_axis_index = axis_names.index("c") if "c" in axis_names else None
    t_axis_index = axis_names.index("t") if "t" in axis_names else None

    if num_channels is not None and channel_axis_index is None:
        raise ValueError(
            f"Input axes {axis_names} contain no 'c' axis; cannot place {num_channels} "
            "output channels."
        )
    if (t_min is not None or t_max is not None) and t_axis_index is None:
        raise ValueError(
            f"Input axes {axis_names} contain no 't' axis; cannot apply --t-min/--t-max."
        )

    # Default: output format matches the input's own storage format, unless
    # explicitly overridden. Uses the shared detect_zarr_format (from
    # ome_zarr_store.py) rather than a locally-rederived heuristic.
    if zarr_format is None:
        zarr_format = detect_zarr_format(src_group)

    datasets = multiscales[0]["datasets"]
    if num_levels is not None:
        if num_levels < 1:
            raise ValueError(f"--num-levels must be >= 1, got {num_levels}")
        if num_levels > len(datasets):
            raise ValueError(
                f"--num-levels {num_levels} exceeds the input's own {len(datasets)} "
                "pyramid level(s)."
            )
        datasets = datasets[:num_levels]
        zattrs["multiscales"][0]["datasets"] = datasets
    level_paths = [ds["path"] for ds in datasets]

    if t_axis_index is not None:
        base_shape_t = src_group[level_paths[0]].shape[t_axis_index]
        t_start = t_min if t_min is not None else 0
        t_end = t_max if t_max is not None else base_shape_t
        if not (0 <= t_start < t_end <= base_shape_t):
            raise ValueError(f"Bad t range [{t_start}, {t_end}) for extent {base_shape_t}.")
        print(f"[Prealloc] Restricting t axis to [{t_start}, {t_end}) of {base_shape_t}")
        # Storing info on time range in this way for now
        zattrs["t_range"] = [t_start, t_end]

    shapes = []
    scales = []
    for ds in datasets:
        src_array = src_group[ds["path"]]
        shape = list(src_array.shape)
        # Only remap the channel axis if --num-channels was explicitly given.
        if num_channels is not None:
            shape[channel_axis_index] = num_channels
        # Only restrict the t axis if --t-min/--t-max was explicitly given.
        if t_axis_index is not None and (t_min is not None or t_max is not None):
            shape[t_axis_index] = t_end - t_start
        shapes.append(tuple(shape))
        scales.append(tuple(ds["coordinateTransformations"][0]["scale"]))

    base_chunks = list(src_group[level_paths[0]].chunks)
    if num_channels is not None:
        base_chunks[channel_axis_index] = num_channels
    if t_axis_index is not None:
        base_chunks[t_axis_index] = min(base_chunks[t_axis_index], shapes[0][t_axis_index])

    builder = ImagePyramidBuilder.from_shapes(
        shapes=shapes,
        base_scale=scales,
        axes=axis_names,
        level_paths=level_paths,
        chunks=tuple(base_chunks),
        data_type=dtype,
        zarr_format=zarr_format,
    )

    dest_root = zarr.open(output_path, mode="w", zarr_format=zarr_format)
    builder.to_zarr(dest_root)

    zattrs["multiscales"][0]["version"] = "0.5" if zarr_format == 3 else "0.4"

    # Needed for cases like ilastik prob stores, which have probs for multiple channels/classes
    if num_channels is not None and "omero" in zattrs:
        zattrs["omero"]["channels"] = _omero_channels(num_channels, channel_labels, dtype)

    dest_root.attrs.update(zattrs)

    for path_name, shape in zip(level_paths, shapes):
        print(f"[Prealloc] Created scale '{path_name}': shape {shape}")
    print(f"[Prealloc] Output Zarr format: v{zarr_format}")


def build_skeleton_explicit(output_path: str, dtype: str,
                             shape: list[int], axis_names: list[str],
                             scale: list[float], chunks: list[int],
                             translation: list[float] | None,
                             zarr_format: int, data_type: str | None,
                             channel_labels: list[str] | None):
    """Mode 2: build a fresh, minimal skeleton from explicit geometry --
    no source zarr to inherit metadata from.

    Design note on what carried over from preallocate_empty_zarr.py / its
    OMEZarrStore.create() and what deliberately didn't:
      - CARRIED OVER: the core idea itself -- build attrs from scratch rather
        than copying a source's zattrs -- is exactly OMEZarrMetadata.create()'s
        approach in ome_zarr_store.py (used by PREALLOCATE_EMPTY_ZARR), and
        the translation convention (`translation = origin_voxels * scale`,
        applied after scale per NGFF's coordinateTransformations ordering)
        matches preallocate_empty_zarr.py's own
        `translation = tuple(origin * scale for origin, scale in zip(global_min, resolution))`.
        That multiplication itself is the CALLER's job (e.g. an updated
        preallocate_empty_zarr.py) -- this function takes the already-computed
        `translation` list directly, so the convention is preserved without
        duplicating that arithmetic here.
      - NOT CARRIED OVER / GENERALIZED BEYOND IT: OMEZarrMetadata.create()
        only ever supports spatial axes -- `all_axis_names = ("z", "y", "x")`,
        sliced to the last `ndim` entries, with no concept of a channel or
        time axis type at all. Since this is meant to be a general community
        tool (not just a labels-store builder), `axis_names`/axis "type"
        here is fully general: any axis name is accepted, and its NGFF
        "type" is inferred (channel/time/space) rather than assumed spatial.
      - NOT CARRIED OVER: chunk size is required explicitly here (--chunks),
        never inferred. preallocate_empty_zarr.py's own rationale for its
        chunk choice (chunk size = tiling stride, so tile writes land on
        chunk boundaries -- see check_grid_alignment) is a caller-side
        decision; this function just takes whatever chunking the caller
        decides is right for its use case.
    """
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

    def _axis_type(name: str) -> str:
        if name == "c":
            return "channel"
        if name == "t":
            return "time"
        return "space"

    # Fresh, minimal multiscales block -- mirrors OMEZarrMetadata.create()'s
    # "build from scratch, don't inherit anything" approach (see docstring
    # above), generalized to any axis set rather than spatial-only.
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


def _omero_channels(num_channels: int, channel_labels: list[str] | None, dtype: str) -> list[dict]:
    """Shared by both modes -- unchanged from the original script's inline logic."""
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
                         help="Mode 1: derive geometry and inherit metadata from this "
                              "existing zarr. Mutually exclusive with --shape/--axes/"
                              "--scale/--chunks (mode 2).")
    parser.add_argument("--output-zarr", required=True)
    parser.add_argument("--dtype", choices=["float32", "uint8", "int32", "uint16", "uint32"],
                         default="float32")
    parser.add_argument("--num-channels", type=int, default=None,
                         help="Mode 1 only: number of output channels. Omit to keep the "
                              "source's own channel count unchanged. In mode 2, channel "
                              "count comes from --shape directly.")
    parser.add_argument("--channel-labels", type=str, default=None,
                         help="Comma-separated channel names. In mode 1, matches "
                              "--num-channels in length. In mode 2, matches --shape's "
                              "channel-axis size.")
    parser.add_argument("--zarr-format", type=int, choices=[2, 3], default=None,
                         help="Output Zarr storage format. Mode 1 defaults to matching "
                              "the input's own format. REQUIRED in mode 2 (no input to "
                              "infer it from).")
    parser.add_argument("--data-type", type=str, default=None,
                         choices=["intensity", "label", "probability"])
    parser.add_argument("--t-min", type=int, default=None,
                         help="Mode 1 only.")
    parser.add_argument("--t-max", type=int, default=None,
                         help="Mode 1 only.")
    parser.add_argument("--num-levels", type=int, default=None,
                         help="Mode 1 only: restrict to the first N pyramid levels. "
                              "Mode 2 is always single-level; passing anything other "
                              "than 1 here in mode 2 is an error.")
    parser.add_argument("--shape", type=str, default=None,
                         help="Mode 2: comma-separated ints, e.g. '10,1024,1024'.")
    parser.add_argument("--axes", type=str, default=None,
                         help="Mode 2: comma-separated axis names matching --shape, "
                              "e.g. 'z,y,x'.")
    parser.add_argument("--scale", type=str, default=None,
                         help="Mode 2: comma-separated per-axis voxel size, matching --shape.")
    parser.add_argument("--chunks", type=str, default=None,
                         help="Mode 2: comma-separated per-axis chunk size, matching --shape. "
                              "Required in mode 2 -- never inferred (e.g. a stitching merge "
                              "target needs this set to the tiling stride so tile writes "
                              "land on chunk boundaries).")
    parser.add_argument("--translation", type=str, default=None,
                         help="Mode 2: comma-separated per-axis physical-space origin, "
                              "matching --shape. Defaults to all zeros if omitted.")

    args = parser.parse_args()

    channel_labels = _parse_csv_strs(args.channel_labels) if args.channel_labels else None
    explicit_mode_args_given = any([args.shape, args.axes, args.scale, args.chunks])

    if args.input_zarr and explicit_mode_args_given:
        parser.error("--input-zarr (mode 1) and --shape/--axes/--scale/--chunks (mode 2) "
                     "are mutually exclusive.")
    if not args.input_zarr and not explicit_mode_args_given:
        parser.error("Provide either --input-zarr (mode 1) or --shape/--axes/--scale/"
                     "--chunks (mode 2).")

    if args.input_zarr:
        if channel_labels and args.num_channels is None:
            parser.error("--channel-labels requires --num-channels to also be set in mode 1.")
        if channel_labels and len(channel_labels) != args.num_channels:
            parser.error(f"--channel-labels has {len(channel_labels)} entries, "
                         f"but --num-channels is {args.num_channels}")
        if args.num_levels is not None and args.num_levels < 1:
            parser.error(f"--num-levels must be >= 1, got {args.num_levels}")

        build_skeleton_from_input(
            args.input_zarr, args.output_zarr, args.dtype,
            args.num_channels, channel_labels, args.zarr_format, args.data_type,
            args.t_min, args.t_max, args.num_levels,
        )
    else:
        if not (args.shape and args.axes and args.scale and args.chunks):
            parser.error("Mode 2 requires --shape, --axes, --scale, and --chunks all set.")
        if args.zarr_format is None:
            parser.error("--zarr-format is required in mode 2 (no input to infer it from).")
        if args.num_levels is not None and args.num_levels != 1:
            parser.error(f"Mode 2 is single-level only; --num-levels must be 1 if given, "
                         f"got {args.num_levels}.")
        if args.num_channels is not None:
            parser.error("--num-channels is not used in mode 2 -- channel count comes "
                         "from --shape directly.")

        shape = _parse_csv_ints(args.shape)
        axis_names = _parse_csv_strs(args.axes)
        scale = _parse_csv_floats(args.scale)
        chunks = _parse_csv_ints(args.chunks)
        translation = _parse_csv_floats(args.translation) if args.translation else None

        build_skeleton_explicit(
            args.output_zarr, args.dtype, shape, axis_names, scale, chunks,
            translation, args.zarr_format, args.data_type, channel_labels,
        )

    print(f"[Prealloc] Done: {args.output_zarr}")


if __name__ == "__main__":
    main()
