#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SAMPLE_INPUT_PREP  } from '../../modules/local/sample_input_prep/main'
include { PYMIF_CONVERSION  } from '../../modules/local/pymif_conversion/main'


workflow {
    samples_to_convert = SAMPLE_INPUT_PREP(
        file(params.exp_dir), 
        params.manifest_name, 
        params.zarr_version
    )
    PYMIF_CONVERSION(samples_to_convert)
}
