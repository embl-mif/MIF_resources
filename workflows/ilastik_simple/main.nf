#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_ZARR_PIPELINE } from '../../subworkflows/local/ilastik_zarr_pipeline/main'

workflow {
    if (!params.input_csv)       error "Missing required parameter: --input_csv"
    if (!params.ilastik_project) error "Missing required parameter: --ilastik_project"

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    ch_ilastik_input = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()

            def input_uri   = Utils.requireColumnValue(row, 'input_image_path')
            def input_image = Utils.resolvePath(input_uri, params.input_csv)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${file(input_uri).simpleName}_${index + 1}"

            def default_zarr_output = params.outdir
                ? "${file(params.outdir)}/probabilities/${dataset_id}_probabilities.zarr"
                : null

            if (!Utils.normalizeOptionalText(row.pixel_class_output) && !default_zarr_output) {
                throw new IllegalArgumentException(
                    "pixel_class_output is blank for dataset_id=${dataset_id} and --outdir is not set. " +
                    "Set an absolute pixel_class_output in the CSV or pass --outdir."
                )
            }

            def zarr_output = Utils.resolveOutput(row, 'pixel_class_output', default_zarr_output, dataset_id)

            def meta = [ id: dataset_id, zarr_output: zarr_output ]

            tuple(meta, input_image, zarr_output)
        }

    ch_project = Channel.value(file(params.ilastik_project))

    ILASTIK_ZARR_PIPELINE(ch_ilastik_input, ch_project)

    ILASTIK_ZARR_PIPELINE.out.probabilities.view()
}
