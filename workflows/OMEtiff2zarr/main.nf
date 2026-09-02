#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SAMPLE_INPUT_PREP  } from '../../modules/local/sample_input_prep/main'
include { PYMIF_CONVERSION  } from '../../modules/local/pymif_conversion/main'


workflow {
    sample_prep_ch = SAMPLE_INPUT_PREP(
        file(params.exp_dir, checkIfExists: true), 
        params.zarr_version,
        params.sample_table
    )
    pymif_inputs_ch = sample_prep_ch.csv
        .splitCsv(header: true)
        .map { row ->
            def ome_tiff = file(row.input)
            def sample_id = file(row.input).simpleName
            def meta = [ id: sample_id, exp_name: file(params.exp_dir).name ]
            return [ meta, row , ome_tiff]
        }
    PYMIF_CONVERSION(pymif_inputs_ch)
}
