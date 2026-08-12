#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRACK_TRACK   } from '../../modules/local/ultrack_track/main'
include { REBUILD_PYRAMID } from '../../modules/local/rebuild_pyramid/main'

workflow {
    def meta = [ id: params.dataset_id ?: "test_dataset" ]

    def prob_zarr = params.probabilities_zarr   // absolute path string
    def raw_image = params.input_zarr

    def ultrack_meta = meta.clone()
    ultrack_meta.foreground_channel = 1   // "cells"
    ultrack_meta.contours_channel   = null
    ultrack_meta.raw_channel        = 0   // brightfield

    ch_ultrack = Channel.of([ ultrack_meta, null, prob_zarr, null, raw_image ])
    ch_config_toml = params.ultrack_config ? file(params.ultrack_config) : []

    ULTRACK_TRACK(ch_ultrack, ch_config_toml)


    REBUILD_PYRAMID(
        ULTRACK_TRACK.out.segments,
        downsample_method: "Sample"
    )
}
