import argparse
import logging

from . import DEFAULT_FEATURES, build_label_path, parse_properties, run_measurement, setup_logging


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract features from OME-Zarr for Nextflow.")
    parser.add_argument("-z", "--input_zarr", required=True, help="Input Zarr path")
    parser.add_argument("-o", "--output_path", required=True, help="Output Parquet path")
    parser.add_argument(
        "--intensity",
        default="0",
        help="Component for intensity (e.g., '0')",
    )
    parser.add_argument(
        "--segmentation",
        default="nuclei",
        help="Segmentation type under the labels group (e.g., 'nuclei', 'membrane')",
    )
    parser.add_argument(
        "--properties",
        default=",".join(DEFAULT_FEATURES),
        help="Comma-separated regionprops properties to extract",
    )
    parser.add_argument(
        "--scheduler",
        choices=["threads", "processes", "single-threaded"],
        default="single-threaded",
        help="Dask scheduler to use for compute; use single-threaded if Nextflow is managing CPU parallelism",
    )
    parser.add_argument(
        "--nextflow-managed",
        action="store_true",
        help="Force single-threaded execution so Nextflow can manage parallelism externally.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()
    setup_logging(args.verbose)

    try:
        features = parse_properties(args.properties)
        label_path = build_label_path(args.segmentation)
        scheduler = "single-threaded" if args.nextflow_managed else args.scheduler
        run_measurement(
            args.input_zarr,
            args.output_path,
            args.intensity,
            label_path,
            features=features,
            scheduler=scheduler,
        )
    except Exception as exc:
        logging.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
