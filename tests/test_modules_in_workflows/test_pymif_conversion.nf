#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PYMIF_CONVERSION } from '../../modules/local/pymif_conversion'

params.test_pymif_csv = "${projectDir}/tests/data/pymif_conversion/input.csv"

workflow {
    input_csv_ch = Channel
        .fromPath(params.test_pymif_csv, checkIfExists: true)
        .map { csv -> [ [ id: csv.baseName, exp_name: 'test_exp' ], csv ] }

    PYMIF_CONVERSION(input_csv_ch)
}
