#!/usr/bin/env python3

import argparse
import copy
from pathlib import Path

import pymif.microscope_manager as mm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a Z maximum-intensity projection using PyMIF."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def remove_spatial_axis(values, index, field):
    if values is None:
        return None

    updated = []

    for level_values in values:
        level_values = tuple(level_values)

        if index >= len(level_values):
            raise ValueError(
                f"Cannot remove spatial axis {index} from "
                f"{field}={level_values}"
            )

        updated.append(
            level_values[:index] + level_values[index + 1 :]
        )

    return updated


def main():
    args = parse_args()

    if args.output.exists():
        raise FileExistsError(f"Output already exists: {args.output}")

    # PyMIF reads all pyramid levels as lazy Dask arrays.
    source = mm.ZarrManager(str(args.input))

    metadata = copy.deepcopy(source.metadata)
    axes = metadata["axes"].lower()

    if "z" not in axes:
        raise ValueError(
            f"Input dataset does not contain a Z axis: axes={axes}"
        )

    z_index = axes.index("z")

    # Apply the MIP independently to every existing pyramid level.
    mip_levels = [
        level.max(axis=z_index)
        for level in source.data
    ]

    output_axes = axes[:z_index] + axes[z_index + 1 :]
    metadata["axes"] = output_axes

    # PyMIF stores scales/units only for spatial axes, rather than all axes.
    spatial_axes = [axis for axis in axes if axis in "zyx"]
    z_spatial_index = spatial_axes.index("z")

    if "scales" in metadata:
        metadata["scales"] = remove_spatial_axis(
            metadata["scales"],
            z_spatial_index,
            "scales",
        )

    if "units" in metadata:
        units = tuple(metadata["units"])

        if len(units) != len(spatial_axes):
            raise ValueError(
                f"Expected {len(spatial_axes)} spatial units for axes "
                f"{axes}, but found {units}"
            )

        metadata["units"] = (
            units[:z_spatial_index]
            + units[z_spatial_index + 1 :]
        )

    # Update level-dependent metadata after removing Z.
    metadata["size"] = [
        tuple(int(size) for size in level.shape)
        for level in mip_levels
    ]

    metadata["chunksize"] = [
        tuple(int(size) for size in level.chunksize)
        for level in mip_levels
    ]

    metadata["dtype"] = str(mip_levels[0].dtype)
    metadata["data_type"] = "intensity"

    # ArrayManager validates the projected arrays against the new axis-aware
    # metadata. PyMIF then writes all OME-NGFF/Zarr metadata.
    output = mm.ArrayManager(mip_levels, metadata)

    output.to_zarr(
        str(args.output),
        overwrite=True,
    )

    source.close()


if __name__ == "__main__":
    main()