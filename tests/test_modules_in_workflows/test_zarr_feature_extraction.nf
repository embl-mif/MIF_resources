#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ZARR_FEATURE_EXTRACTION } from '../../modules/local/zarr_feature_extraction'

def test_zarr = params.test_zarr
println "TEST_ZARR=${test_zarr}"

params.intensity = "0"
params.labels    = "nuclei"

workflow {
    println "DEBUG: outdir    = ${params.outdir}"
    println "DEBUG: intensity = ${params.intensity}"
    println "DEBUG: labels    = ${params.labels}"

    zarr_ch = Channel
        .fromPath(test_zarr, type: 'dir', checkIfExists: true)
        .map { it -> [ [ id: it.baseName, exp_name: 'test_exp' ], it ] }

    ZARR_FEATURE_EXTRACTION (
        zarr_ch,
        params.intensity,
        params.labels
    )
}
