include { PREALLOCATE_ZARR } from '../../../modules/local/preallocate_zarr/main'
include { RUN_ULTRACK      } from '../../../modules/local/run_ultrack/main'
include { POPULATE_PYRAMID } from '../../../modules/local/populate_pyramid/main'

workflow ULTRACK_ZARR_PIPELINE {

    take:
    ch_ultrack_input    // channel: [ meta, prob_zarr, raw_image ]
                        //   meta carries foreground_channel/raw_channel/t_min/t_max/
                        //   min_area/max_area/max_distance/solution_gap/time_limit etc.
    val_config_toml     // Channel.value(path(config.toml)) or Channel.value([])

    main:
    ch_versions = Channel.empty()
    def repoRootDir = workflow.projectDir.getParent().getParent()
    def noTilesCsv     = file("${repoRootDir}/assets/NO_FILE")
    def noImageSentinel = file("${repoRootDir}/assets/NO_FILE")
    // This pipeline never has labels_zarr/contours_zarr of its own -- both
    // are passed to RUN_ULTRACK as this same sentinel, same convention as
    // CELLPOSE_SAM_ZARR_MIF's optional model input.
    def noOptionalZarr  = file("${repoRootDir}/assets/NO_FILE")

    // raw_image and prob_zarr can each be a local file or a remote URI --
    // split once here so every downstream process gets the path(local)/
    // val(uri) pair it needs for correct container bind-mount visibility.
    ch_split_ultrack_input = ch_ultrack_input
        .map { meta, prob_zarr, raw_image ->
            def (localRawImage, rawImageUri) = Utils.splitImageInput(raw_image, noImageSentinel)
            def (localProbZarr, probZarrUri) = Utils.splitImageInput(prob_zarr, noImageSentinel)
            tuple(meta.id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri)
        }

    ch_prealloc_input = ch_split_ultrack_input
        .map { id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri ->
            def prealloc_meta = meta.clone()
            prealloc_meta.num_channels = 1
            prealloc_meta.dtype        = "int32"
            prealloc_meta.data_type    = "label"
            tuple(prealloc_meta, localRawImage, rawImageUri, noTilesCsv, meta.output_zarr)
        }

    PREALLOCATE_ZARR(ch_prealloc_input)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_ultrack = ch_split_ultrack_input
        .combine(
            PREALLOCATE_ZARR.out.zarr.map { meta, zarr -> tuple(meta.id, zarr) },
            by: 0
        )
        .map { id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri, prealloc_zarr ->
            // labels_zarr/contours_zarr: not used by this pipeline -- sentinel,
            // no uri. foreground_zarr: the probability map. raw_zarr: the
            // original image (used only to derive contours if contours_zarr
            // is absent, per RUN_ULTRACK's own logic). output_tracked_zarr:
            // the store PREALLOCATE_ZARR just built, staged in as `path` so
            // RUN_ULTRACK's container can see it and -resume can verify it.
            tuple(
                meta,
                noOptionalZarr, null,
                localProbZarr,  probZarrUri,
                noOptionalZarr, null,
                localRawImage,  rawImageUri,
                prealloc_zarr
            )
        }

    RUN_ULTRACK(ch_ultrack, val_config_toml)
    ch_versions = ch_versions.mix(RUN_ULTRACK.out.versions)

    // downsample method = "Sample" for integer labels, not "Average" --
    // averaging would blend adjacent cell IDs into meaningless values.
    ch_ready_for_pyramid = RUN_ULTRACK.out.tracked_zarr
        .map { meta, output_zarr -> tuple(meta, output_zarr, "Sample") }

    POPULATE_PYRAMID(ch_ready_for_pyramid)
    ch_versions = ch_versions.mix(POPULATE_PYRAMID.out.versions)

    emit:
    tracked_zarr = POPULATE_PYRAMID.out.zarr
    versions     = ch_versions
}
