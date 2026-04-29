import logging
import os
import sys

import dask
import dask.array as da
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from skimage.measure import regionprops_table

DEFAULT_FEATURES = ["label", "area", "centroid", "intensity_mean", "intensity_max"]


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(levelname)s] %(message)s",
        stream=sys.stdout,
    )


def parse_properties(properties: str) -> list[str]:
    return [item.strip() for item in properties.split(",") if item.strip()]


def validate_zarr_components(input_zarr: str, intensity_path: str, label_path: str) -> zarr.Group:
    if not os.path.exists(input_zarr):
        raise FileNotFoundError(f"Input Zarr path does not exist: {input_zarr}")

    try:
        group = zarr.open_group(input_zarr, mode="r")
    except Exception as exc:
        raise RuntimeError(f"Could not open Zarr group at {input_zarr}: {exc}") from exc

    missing = []
    for component in (intensity_path, label_path):
        try:
            group[component]
        except KeyError:
            missing.append(component)

    if missing:
        raise ValueError(f"Missing Zarr components: {', '.join(missing)}")

    return group


def make_output_parent(output_path: str) -> None:
    parent_dir = os.path.dirname(os.path.abspath(output_path))
    if parent_dir and not os.path.isdir(parent_dir):
        os.makedirs(parent_dir, exist_ok=True)


def build_label_path(segmentation_type: str) -> str:
    if segmentation_type.startswith("labels/") or segmentation_type.startswith("labels\\"):
        return segmentation_type
    return f"labels/{segmentation_type}"


def prepare_intensity_image(intensity_t: np.ndarray, labels_shape: tuple[int, ...]) -> np.ndarray:
    intensity = np.asarray(intensity_t)

    if intensity.ndim > len(labels_shape):
        intensity = np.moveaxis(intensity, 0, -1)

    if intensity.shape[: len(labels_shape)] != labels_shape:
        raise ValueError(
            "Intensity image shape %s is incompatible with labels shape %s"
            % (intensity.shape, labels_shape)
        )

    return intensity


def run_measurement(
    input_zarr: str,
    output_parquet: str,
    intensity_path: str,
    label_path: str,
    features: list[str] | None = None,
    scheduler: str = "single-threaded",
) -> str:
    if features is None:
        features = DEFAULT_FEATURES

    dask.config.set(scheduler=scheduler)
    logger = logging.getLogger(__name__)

    validate_zarr_components(input_zarr, intensity_path, label_path)
    make_output_parent(output_parquet)

    intensity_stack = da.from_zarr(input_zarr, component=intensity_path)
    label_stack = da.from_zarr(input_zarr, component=label_path)

    if label_stack.shape[0] != intensity_stack.shape[0]:
        raise ValueError(
            "Zarr intensity and label stacks must have the same number of timepoints. "
            f"Found {intensity_stack.shape[0]} and {label_stack.shape[0]}"
        )

    n_tps = label_stack.shape[0]
    logger.info("Processing %d timepoint(s)", n_tps)

    writer = None
    try:
        for t in range(n_tps):
            logger.info("Processing timepoint %d/%d", t + 1, n_tps)
            labels_t = label_stack[t].compute()
            intensity_t = intensity_stack[t].compute()

            intensity_t = prepare_intensity_image(intensity_t, labels_t.shape)
            props = regionprops_table(labels_t, intensity_t, properties=features)

            df = pd.DataFrame(props)
            df["t"] = t

            table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_parquet, table.schema, compression="snappy")

            if len(df):
                writer.write_table(table)
            else:
                logger.debug("No objects detected for timepoint %d", t)

        if writer is None:
            empty_df = pd.DataFrame(columns=features + ["t"])
            pq.write_table(pa.Table.from_pandas(empty_df, preserve_index=False), output_parquet, compression="snappy")
    finally:
        if writer is not None:
            writer.close()

    logger.info("Successfully saved to %s", output_parquet)
    return output_parquet
