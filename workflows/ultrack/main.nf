#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREALLOCATE_ZARR } from '../../modules/local/preallocate_zarr/main'
include { RUN_ULTRACK          } from '../../modules/local/run_ultrack/main'
include { POPULATE_PYRAMID } from '../../modules/local/populate_pyramid/main'

workflow {
    def pubdir = "${file(params.outdir)}/tracking"
    def meta =  [ id: "test_dataset" ]

    def prob_zarr = params.probabilities_zarr   // absolute path string
    def raw_image = params.input_zarr
    def output_tracked_zarr = "${pubdir}/${meta.id}_tracked.zarr"

    def prealloc_meta = meta.clone()
    prealloc_meta.num_channels = 1
    prealloc_meta.dtype        = "int32"
    prealloc_meta.data_type    = "label"

    ch_prealloc_input = Channel.of([ prealloc_meta, raw_image, output_tracked_zarr ])
    PREALLOCATE_ZARR(ch_prealloc_input)

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

    ch_ultrack = PREALLOCATE_ZARR.out.empty_zarr
        .map { m, prealloc_zarr -> tuple(ultrack_meta, null, prob_zarr, null, raw_image, prealloc_zarr) }
    ch_config_toml = params.ultrack_config ? file(params.ultrack_config) : []

    RUN_ULTRACK(ch_ultrack, ch_config_toml)
    ch_ready_for_pyramid = RUN_ULTRACK.out.tracked_zarr
        .map { m, output_zarr -> tuple(m, output_zarr, "Sample") } // downsample methos = sample for int

    POPULATE_PYRAMID(ch_ready_for_pyramid)
}
