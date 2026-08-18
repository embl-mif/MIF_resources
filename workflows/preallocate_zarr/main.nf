#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREALLOCATE_ZARR } from '../../modules/local/preallocate_zarr/main'

workflow {
    def meta = [
        id: "test_dataset",
        num_channels: 3,
        channel_labels: ["background", "cells", "crystals"],
        zarr_format: 3
    ]

    def input_zarr  = params.input_zarr
    def output_zarr = "${params.outdir}/${meta.id}_test.zarr"

    ch_prealloc_input = Channel.of([ meta, input_zarr, output_zarr ])

    PREALLOCATE_ZARR(ch_prealloc_input)

    PREALLOCATE_ZARR.out.empty_zarr.view()
}
