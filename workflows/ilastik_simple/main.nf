#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILASTIK_ZARR_PIPELINE } from '../../subworkflows/local/ilastik_zarr_pipeline/main'

workflow {
    if (!params.input_zarr)      error "Missing required parameter: --input_zarr"
    if (!params.ilastik_project) error "Missing required parameter: --ilastik_project"
    if (!params.outdir)          error "Missing required parameter: --outdir"

    def meta = [ id: params.dataset_id ?: "test_dataset" ]
    def output_zarr = "${params.outdir}/probabilities/${meta.id}_probabilities.zarr"

    ch_ilastik_input = Channel.of([ meta, params.input_zarr, output_zarr ])
    ch_project = Channel.value(file(params.ilastik_project))

    ILASTIK_ZARR_PIPELINE(ch_ilastik_input, ch_project)

    ILASTIK_ZARR_PIPELINE.out.probabilities.view()
}
