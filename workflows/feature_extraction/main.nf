#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_FEATURE_EXTRACTION } from '../../modules/local/pymif_feature_extraction/main'

workflow {

    def default_features_dir = "${file(params.default_outdir)}/measurements"

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    feature_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def intensity_zarr = file(Utils.requireColumnValue(row, 'intensity'), checkIfExists: true)
            def labels_zarr = file(Utils.requireColumnValue(row, 'labels'), checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${intensity_zarr.simpleName}_${index + 1}"

            def requested_features_output = Utils.normalizeOptionalText(row.features_output)
            if (requested_features_output && !file(requested_features_output).isAbsolute()) {
                error(
                    "features_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_features_output}'. Leave the column blank to use the " +
                    "default location under default_outdir instead."
                )
            }
            def resolved_features_output = requested_features_output
                ?: "${default_features_dir}/${dataset_id}_features.csv"

            def intensity_resolution_level = Utils.normalizeOptionalText(row.intensity_resolution_level)
            def labels_resolution_level = Utils.normalizeOptionalText(row.labels_resolution_level)

            def meta = row + [
                id                                       : dataset_id,
                features_output_name                     : resolved_features_output,
                features_intensity_resolution_level       : intensity_resolution_level
                    ? intensity_resolution_level as int
                    : (params.intensity_resolution_level != null ? params.intensity_resolution_level as int : 0),
                features_labels_resolution_level          : labels_resolution_level
                    ? labels_resolution_level as int
                    : (params.labels_resolution_level != null ? params.labels_resolution_level as int : 0),
                features_channels                         : Utils.normalizeOptionalText(row.channels),
                features_extra_properties                 : Utils.normalizeOptionalText(row.extra_properties),
                features_extra_intensity_properties       : Utils.normalizeOptionalText(row.extra_intensity_properties),
            ]

            tuple(meta, intensity_zarr, labels_zarr)
        }

    PYMIF_FEATURE_EXTRACTION(feature_inputs_ch)

}
