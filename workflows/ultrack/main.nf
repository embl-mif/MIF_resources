#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRACK_ZARR_PIPELINE } from '../../subworkflows/local/ultrack_zarr_pipeline/main'

workflow {
    if (!params.input_csv) error "Missing required parameter: --input_csv"

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    ch_ultrack_input = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()

            def prob_uri  = Utils.requireColumnValue(row, 'probabilities_zarr')
            def input_uri = Utils.requireColumnValue(row, 'input_zarr')

            def probabilities_zarr = Utils.resolvePath(prob_uri, params.input_csv)
            def input_zarr         = Utils.resolvePath(input_uri, params.input_csv)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${file(input_uri).simpleName}_${index + 1}"

            def default_zarr_output = params.outdir
                ? "${file(params.outdir)}/tracking/${dataset_id}_tracked.zarr"
                : null

            if (!Utils.normalizeOptionalText(row.ultrack_output) && !default_zarr_output) {
                throw new IllegalArgumentException(
                    "ultrack_output is blank for dataset_id=${dataset_id} and --outdir is not set. " +
                    "Set an absolute ultrack_output in the CSV or pass --outdir."
                )
            }

            def zarr_output = Utils.resolveOutput(row, 'ultrack_output', default_zarr_output, dataset_id)

            tuple(row + [id: dataset_id, zarr_output: zarr_output], probabilities_zarr, input_zarr)
        }

    val_config_toml = params.ultrack_config ? Channel.value(file(params.ultrack_config)) : Channel.value([])

    ULTRACK_ZARR_PIPELINE(ch_ultrack_input, val_config_toml)

    ULTRACK_ZARR_PIPELINE.out.tracked_zarr.view()
}
