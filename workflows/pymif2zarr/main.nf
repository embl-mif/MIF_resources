#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PYMIF_CONVERSION } from '../../modules/local/pymif_conversion/main'


workflow {
    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    def resolveOutput = { value, fallback ->
        def v = Utils.normalizeOptionalText(value)
        !v ? fallback : (v.startsWith('/') ? v : file(v).name)
    }

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def resolved_zarr_output = resolveOutput(row.output, "${dataset_id}.ome.zarr")
            
            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)
}
