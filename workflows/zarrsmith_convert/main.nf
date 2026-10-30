#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ZARRSMITH_CONVERSION } from '../../modules/local/zarrsmith_conversion/main'


workflow {
    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    zarrsmith_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def resolved_zarr_output = Utils.resolveOutput(row, 'output', "${dataset_id}.ome.zarr", dataset_id)

            def meta = row + [
                id                       : dataset_id,
                zarr_output         : resolved_zarr_output,
            ]

            tuple(meta, input_dataset)
        }

    ZARRSMITH_CONVERSION(zarrsmith_inputs_ch)
}
