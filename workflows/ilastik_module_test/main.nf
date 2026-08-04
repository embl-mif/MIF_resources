#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ILASTIK_SUBREGION_PIXEL_CLASS  } from '../../modules/local/ilastik_subregion_pixel_class/main'
include { REBUILD_PYRAMID               } from '../../modules/local/rebuild_pyramid/main'

workflow {
    def pubdir = "${file(params.outdir)}/probabilities"
    
    def meta = [ id: params.dataset_id ?: "test_dataset" ]
    def raw_zarr = file(params.input_zarr)
    def master_output_zarr = "${pubdir}/${meta.id}_probabilities.zarr"

    // Emit single tuple directly
    ilastik_input_ch = Channel.of([ meta, raw_zarr, master_output_zarr ])
    ch_project       = Channel.value(file(params.ilastik_project))

    ILASTIK_SUBREGION_PIXEL_CLASS(ilastik_input_ch, ch_project)

    ch_ready_for_pyramid = ILASTIK_SUBREGION_PIXEL_CLASS.out.zarr
        .collect(flat: false)
        .map { it[0] }

    REBUILD_PYRAMID(ch_ready_for_pyramid)
}
