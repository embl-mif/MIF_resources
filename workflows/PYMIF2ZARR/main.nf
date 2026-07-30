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

            def output_dataset = file(row.output, checkIfExists: false)
            if( !output_dataset.isAbsolute() ) {
                error "CSV column 'output' must be an absolute path, got: ${row.output}"
            }

            def meta = [
                id: input_dataset.simpleName,
                outdir: output_dataset.parent ?: input_dataset.parent
            ]

            def row2 = row + [
                output_name: output_dataset.name,
            ]
            return [meta, row2, input_dataset]
        }

    PYMIF_CONVERSION(pymif_inputs_ch)
}
