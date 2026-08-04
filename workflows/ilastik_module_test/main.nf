#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ILASTIK_ZARR  } from '../../modules/local/ilastik_zarr/main'
include { REBUILD_PYRAMID               } from '../../modules/local/rebuild_pyramid/main'

workflow {
    def pubdir = "${file(params.outdir)}/probabilities"
    
    def meta = [ id: params.dataset_id ?: "test_dataset" ]
    def raw_zarr = file(params.input_zarr)
    def master_output_zarr = "${pubdir}/${meta.id}_probabilities.zarr"

    // Emit single tuple directly
    ilastik_input_ch = Channel.of([ meta, raw_zarr, master_output_zarr ])
    ch_project       = Channel.value(file(params.ilastik_project))

    ILASTIK_ZARR(ilastik_input_ch, ch_project)

    ch_ready_for_pyramid = ILASTIK_ZARR.out.zarr
        .map { tileMeta, output_zarr -> tuple(tileMeta.id, tileMeta, output_zarr) }
        .groupTuple(by: 0)
        .map { id, metas, zarrs -> tuple(metas[0], zarrs[0]) }

    REBUILD_PYRAMID(ch_ready_for_pyramid)
}
