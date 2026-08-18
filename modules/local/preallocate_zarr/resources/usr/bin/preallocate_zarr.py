#!/usr/bin/env python3
"""Preallocate an empty OME-Zarr skeleton from a source image's metadata,
optionally remapping the channel axis to a different count/labels (e.g.
ilastik probability outputs). Uses ngio's ImagePyramidBuilder for v2/v3
dual support. See https://github.com/BioVisionCenter/ngio
"""

import argparse

import zarr
from ngio.common._pyramid import ImagePyramidBuilder


def _detect_input_zarr_format(root) -> int:
    """v3 stores expose zarr_format via root.metadata; v2 stores don't."""
    metadata = getattr(root, "metadata", None)
    if metadata is None:
        return 2
    return getattr(metadata, "zarr_format", 2)


def build_skeleton(input_path: str, output_path: str, dtype: str,
                    num_channels: int | None, channel_labels: list[str] | None,
                    zarr_format: int | None):
    src_group = zarr.open(input_path, mode="r")
    zattrs = dict(src_group.attrs)
    multiscales = zattrs.get("multiscales", [])
    if not multiscales:
        raise ValueError(f"Input Zarr '{input_path}' lacks 'multiscales' metadata.")

    axes_meta = multiscales[0].get("axes", [])
    axis_names = tuple(
        (axis.get("name") if isinstance(axis, dict) else axis) for axis in axes_meta
    )
    channel_axis_index = axis_names.index("c") if "c" in axis_names else None

    if num_channels is not None and channel_axis_index is None:
        raise ValueError(
            f"Input axes {axis_names} contain no 'c' axis; cannot place {num_channels} "
            "output channels."
        )

    # Default: output format matches the input's own storage format, unless
    # explicitly overridden.
    if zarr_format is None:
        zarr_format = _detect_input_zarr_format(src_group)

    datasets = multiscales[0]["datasets"]
    level_paths = [ds["path"] for ds in datasets]

    shapes = []
    scales = []
    for ds in datasets:
        src_array = src_group[ds["path"]]
        shape = list(src_array.shape)
        # Only remap the channel axis if --num-channels was explicitly given.
        if num_channels is not None:
            shape[channel_axis_index] = num_channels
        shapes.append(tuple(shape))
        scales.append(tuple(ds["coordinateTransformations"][0]["scale"]))

    base_chunks = list(src_group[level_paths[0]].chunks)
    if num_channels is not None:
        base_chunks[channel_axis_index] = num_channels

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
        max_val = 255.0 if dtype == "uint8" else 1.0
        palette = ["FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "FFFFFF"]
        labels = channel_labels or [f"channel_{i}" for i in range(num_channels)]
        zattrs["omero"]["channels"] = [
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

    dest_root.attrs.update(zattrs)

    for path_name, shape in zip(level_paths, shapes):
        print(f"[Prealloc] Created scale '{path_name}': shape {shape}")
    print(f"[Prealloc] Output Zarr format: v{zarr_format}")


def main():
    parser = argparse.ArgumentParser(
        description="Preallocate an empty OME-Zarr skeleton, optionally remapping the channel axis."
    )
    parser.add_argument("--input-zarr", required=True)
    parser.add_argument("--output-zarr", required=True)
    parser.add_argument("--dtype", choices=["float32", "uint8", "int32", "uint16", "uint32"],
                         default="float32")
    parser.add_argument("--num-channels", type=int, default=None,
                         help="Number of output channels. Omit to keep the source's own "
                              "channel count unchanged (no remapping).")
    parser.add_argument("--channel-labels", type=str, default=None,
                         help="Comma-separated channel names, matching --num-channels in length. "
                              "Only used if --num-channels is given.")
    parser.add_argument("--zarr-format", type=int, choices=[2, 3], default=None,
                         help="Output Zarr storage format. Defaults to matching the input's own "
                              "format. Set explicitly to force v2 or v3 output regardless of input.")

    args = parser.parse_args()
    channel_labels = args.channel_labels.split(",") if args.channel_labels else None
    if channel_labels and args.num_channels is None:
        parser.error("--channel-labels requires --num-channels to also be set.")
    if channel_labels and len(channel_labels) != args.num_channels:
        parser.error(f"--channel-labels has {len(channel_labels)} entries, "
                     f"but --num-channels is {args.num_channels}")

    build_skeleton(args.input_zarr, args.output_zarr, args.dtype,
                    args.num_channels, channel_labels, args.zarr_format)
    print(f"[Prealloc] Done: {args.output_zarr}")


if __name__ == "__main__":
    main()
