#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ZARR_FEATURE_EXTRACTION } from '../../modules/local/zarr_feature_extraction/main.nf'

def test_zarr = '/home/vivanco/mif_resources/tests/data/microscopy/A1_F1.zarr'
println "TEST_ZARR=${test_zarr}"

workflow {

    zarr_ch = Channel
        .fromPath(test_zarr)
    zarr_ch.view()
}
