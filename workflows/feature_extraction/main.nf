#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ZARR_FEATURE_EXTRACTION } from '../../modules/local/zarr_feature_extraction'


workflow {
    println "DEBUG: outdir    = ${params.outdir}"
    println "DEBUG: intensity = ${params.intensity}"
    println "DEBUG: labels    = ${params.labels}"

    zarr_ch = Channel
        .fromPath(params.test_zarr, type: 'dir', checkIfExists: true)
        .map { it -> [ [ id: it.baseName, exp_name: 'test_exp' ], it ] }

    ZARR_FEATURE_EXTRACTION (
        zarr_ch,
        params.intensity,
        params.labels
    )
}

