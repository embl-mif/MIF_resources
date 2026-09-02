#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PYMIF_CONVERSION } from '../../modules/local/pymif_conversion/main'


workflow {
    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = row.dataset_id?.toString()?.trim()
                ?: input_dataset.simpleName
            
            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : row.output,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)
}
