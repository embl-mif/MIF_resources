include { EXTRACT_ILP_CLASSES   } from '../../../modules/local/extract_ilp_classes/main'
include { PREALLOCATE_ZARR      } from '../../../modules/local/preallocate_zarr/main'
include { ILASTIK_ZARR          } from '../../../modules/local/ilastik_zarr/main'
include { POPULATE_PYRAMID      } from '../../../modules/local/populate_pyramid/main'

workflow ILASTIK_ZARR_PIPELINE {

    take:
    ilastik_input_ch   // channel: [ meta, raw_zarr, output_zarr ]
                        //   -- one entry per tile (tiled case)
                        //   -- or one entry per whole dataset (non-tiled case)
    ch_project          // Channel.value(path(ilp))

    main:
    ch_versions = Channel.empty()

    EXTRACT_ILP_CLASSES(ch_project)
    ch_versions = ch_versions.mix(EXTRACT_ILP_CLASSES.out.versions)

    ch_channel_info = EXTRACT_ILP_CLASSES.out.classes
        .map { json -> new groovy.json.JsonSlurper().parseText(json.text) }

    ch_prealloc_input = ilastik_input_ch
        .combine(ch_channel_info)
        .map { meta, raw_zarr, output_zarr, channelInfo ->
            def prealloc_meta = meta.clone()
            prealloc_meta.num_channels   = channelInfo.num_channels
            prealloc_meta.channel_labels = channelInfo.channel_labels
            prealloc_meta.dtype          = "float32"
            tuple(prealloc_meta, raw_zarr, output_zarr)
        }

    ch_prealloc_unique = ch_prealloc_input
        .map { meta, raw_zarr, output_zarr -> tuple(output_zarr, meta, raw_zarr) }
        .unique { it[0] }
        .map { output_zarr, meta, raw_zarr -> tuple(meta, raw_zarr, output_zarr) }

    PREALLOCATE_ZARR(ch_prealloc_unique)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_ilastik_ready = ilastik_input_ch
        .map { meta, raw_zarr, output_zarr -> tuple(output_zarr, meta, raw_zarr) }
        .combine(
            PREALLOCATE_ZARR.out.empty_zarr.map { m, zarr -> tuple(zarr, true) },
            by: 0
        )
        .map { output_zarr, meta, raw_zarr, ready -> tuple(meta, raw_zarr, output_zarr) }

    ILASTIK_ZARR(ch_ilastik_ready, ch_project)
    ch_versions = ch_versions.mix(ILASTIK_ZARR.out.versions)

    // Group per dataset (meta.id) so pyramid rebuild waits for every tile
    ch_ready_for_pyramid = ILASTIK_ZARR.out.zarr
        .map { meta, output_zarr -> tuple(meta.id, meta, output_zarr) }
        .groupTuple(by: 0)
        .map { id, metas, zarrs -> tuple(metas[0], zarrs[0], "Average") }

    POPULATE_PYRAMID(ch_ready_for_pyramid)
    ch_versions = ch_versions.mix(POPULATE_PYRAMID.out.versions)

    emit:
    probabilities = POPULATE_PYRAMID.out.zarr
    versions      = ch_versions
}
