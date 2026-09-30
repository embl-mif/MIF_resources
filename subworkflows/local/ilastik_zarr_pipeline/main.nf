include { EXTRACT_ILASTIK_CLASSES } from '../../../modules/local/extract_ilastik_classes/main'
include { PREALLOCATE_ZARR        } from '../../../modules/local/preallocate_zarr/main'
include { ILASTIK_ZARR            } from '../../../modules/local/ilastik_zarr/main'
include { POPULATE_PYRAMID        } from '../../../modules/local/populate_pyramid/main'

workflow ILASTIK_ZARR_PIPELINE {

    take:
    ilastik_input_ch   // channel: [ meta, raw_zarr, output_zarr ]
                        //   -- one entry per tile (tiled case)
                        //   -- or one entry per whole dataset (non-tiled case)
    ch_project          // Channel.value(path(ilp))

    main:
    ch_versions = Channel.empty()

    def repoRootDir = workflow.projectDir.getParent().getParent()
    def noTilesCsv      = file("${repoRootDir}/assets/NO_FILE")
    def noImageSentinel = file("${repoRootDir}/assets/NO_FILE")

    EXTRACT_ILASTIK_CLASSES(ch_project)
    ch_versions = ch_versions.mix(EXTRACT_ILASTIK_CLASSES.out.versions)

    ch_channel_info = EXTRACT_ILASTIK_CLASSES.out.classes
        .map { json -> new groovy.json.JsonSlurper().parseText(json.text) }

    // raw_zarr can be a local file or a remote URI -- split once so every
    // downstream process gets the path(local)/val(uri) pair it needs.
    ch_split_ilastik_input = ilastik_input_ch
        .map { meta, raw_zarr, output_zarr ->
            def (localRawZarr, rawZarrUri) = Utils.splitImageInput(raw_zarr, noImageSentinel)
            tuple(meta, localRawZarr, rawZarrUri, output_zarr)
        }

    ch_prealloc_input = ch_split_ilastik_input
        .combine(ch_channel_info)
        .map { meta, localRawZarr, rawZarrUri, output_zarr, channelInfo ->
            def prealloc_meta = meta.clone()
            prealloc_meta.num_channels   = channelInfo.num_channels
            prealloc_meta.channel_labels = channelInfo.channel_labels
            prealloc_meta.dtype          = "float32"
            prealloc_meta.t_min = null
            prealloc_meta.t_max = null
            prealloc_meta.z_min = null
            prealloc_meta.z_max = null
            prealloc_meta.y_min = null
            prealloc_meta.y_max = null
            prealloc_meta.x_min = null
            prealloc_meta.x_max = null
            prealloc_meta.output_zarr_key = output_zarr
            tuple(prealloc_meta, localRawZarr, rawZarrUri, output_zarr)
        }

    ch_prealloc_unique = ch_prealloc_input
        .map { meta, localRawZarr, rawZarrUri, output_zarr -> tuple(output_zarr, meta, localRawZarr, rawZarrUri) }
        .unique { it[0] }
        .map { output_zarr, meta, localRawZarr, rawZarrUri -> tuple(meta, localRawZarr, rawZarrUri, noTilesCsv, output_zarr) }

    PREALLOCATE_ZARR(ch_prealloc_unique)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_ilastik_ready = ch_split_ilastik_input
        .map { meta, localRawZarr, rawZarrUri, output_zarr -> tuple(output_zarr, meta, localRawZarr, rawZarrUri) }
        .combine(
            PREALLOCATE_ZARR.out.zarr.map { m, zarr -> tuple(m.output_zarr_key, zarr) },
            by: 0
        )
        .map { output_zarr, meta, localRawZarr, rawZarrUri, prealloc_zarr ->
            tuple(meta, localRawZarr, rawZarrUri, prealloc_zarr)
        }

    ILASTIK_ZARR(ch_ilastik_ready, ch_project)
    ch_versions = ch_versions.mix(ILASTIK_ZARR.out.versions)

    // Group per dataset (meta.id) so pyramid rebuild waits for every tile.
    // Probabilities are continuous float values, not label IDs, so "Average"
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
