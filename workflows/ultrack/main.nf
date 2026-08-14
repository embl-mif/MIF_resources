#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRACK   } from '../../modules/local/ultrack/main'
include { REBUILD_PYRAMID } from '../../modules/local/rebuild_pyramid/main'

workflow {
    def pubdir = "${file(params.outdir)}/tracking"
    def meta =  [ id: "test_dataset" ]

    def prob_zarr = params.probabilities_zarr   // absolute path string
    def raw_image = params.input_zarr
    def output_tracked_zarr = "${pubdir}/${meta.id}_tracked.zarr"

    def ultrack_meta = meta.clone()
    ultrack_meta.foreground_channel = 1   // "cells"
    ultrack_meta.contours_channel   = null
    ultrack_meta.raw_channel        = 0   // brightfield
    ultrack_meta.t_min          = 0
    ultrack_meta.t_max          = 5
    ultrack_meta.min_area       = 200
    ultrack_meta.max_area       = 20000
    ultrack_meta.max_distance   = 10.0
    ultrack_meta.solution_gap   = 0.1
    ultrack_meta.time_limit     = 600
    ch_ultrack = Channel.of([ ultrack_meta, null, prob_zarr, null, raw_image ])
    ch_config_toml = params.ultrack_config ? file(params.ultrack_config) : []

    ULTRACK(ch_ultrack, ch_config_toml)
    ch_ready_for_pyramid = ULTRACK.out.segments
        .map { m, output_zarr -> tuple(m, output_zarr, "Sample") } // downsample methos = sample for int

    REBUILD_PYRAMID(ch_ready_for_pyramid)
}
