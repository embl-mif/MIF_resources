#!/usr/bin/env nextflow

//import module
include { FEATURE_EXTRACTION } from '../../modules/local/pymif/zarr_feature_extraction/main.nf'

workflow {

    def meta = [exp_name: 'test_exp', well: 'A1', field: 'F1']

    def zarr_path = file(params.test_zarr)

    ch_input = Channel.of([ zarr_path, meta ])

    FEATURE_EXTRACTION(ch_input)
}
