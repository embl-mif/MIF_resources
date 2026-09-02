#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRACK_ZARR_PIPELINE } from '../../subworkflows/local/ultrack_zarr_pipeline/main'

workflow {
    if (!params.probabilities_zarr) error "Missing required parameter: --probabilities_zarr"
    if (!params.input_zarr)         error "Missing required parameter: --input_zarr"
    if (!params.outdir)             error "Missing required parameter: --outdir"

    def meta = [
        id                 : params.dataset_id ?: "test_dataset",
        output_dir         : "${file(params.outdir)}/tracking",
        foreground_channel : 1,      // "cells"
        contours_channel   : null,
        raw_channel        : 0,      // brightfield
        t_min              : 0,
        t_max              : 5,
        min_area           : 200,
        max_area            : 20000,
        max_distance        : 10.0,
        solution_gap         : 0.1,
        time_limit           : 600,
    ]
    ch_ultrack_input = Channel.of([ meta, params.probabilities_zarr, params.input_zarr ])
    val_config_toml   = params.ultrack_config ? Channel.value(file(params.ultrack_config)) : Channel.value([])

    ULTRACK_ZARR_PIPELINE(ch_ultrack_input, val_config_toml)

    ULTRACK_ZARR_PIPELINE.out.tracked_zarr.view()
}
