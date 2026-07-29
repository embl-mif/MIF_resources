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

            def output_path = new File(row.output?.toString() ?: '')
            if( !output_path.isAbsolute() ) {
                error "CSV column 'output' must be an absolute path, got: ${row.output}"
            }

            def meta = [
                id: input_dataset.simpleName,
                sample_name: row.sample_name ?: input_dataset.simpleName
            ]

            def row2 = row + [
                output_name: output_path.name,
                output_dir : output_path.parent
            ]

            return [meta, row2, input_dataset]
        }

    PYMIF_CONVERSION(pymif_inputs_ch)
}
