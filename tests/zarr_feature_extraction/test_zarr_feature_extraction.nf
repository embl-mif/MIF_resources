#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ZARR_FEATURE_EXTRACTION } from '../../modules/local/pymif/zarr_feature_extraction/main.nf'

workflow {

    zarr_ch = Channel
        .fromPath(params.test_zarr)
        .map { file -> 
            tuple([id: file.baseName], file) 
        }

    ZARR_FEATURE_EXTRACTION(
        zarr_ch, 
        "0", 
        "labels/nuclei"
    )
}
