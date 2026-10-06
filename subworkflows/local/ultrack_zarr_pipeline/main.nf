include { PREALLOCATE_ZARR } from '../../../modules/local/preallocate_zarr/main'
include { RUN_ULTRACK      } from '../../../modules/local/run_ultrack/main'
include { POPULATE_PYRAMID } from '../../../modules/local/populate_pyramid/main'

workflow ULTRACK_ZARR_PIPELINE {

    take:
    ch_ultrack_input
    val_config_toml

    main:
    ch_versions = Channel.empty()
    def repoRootDir = workflow.projectDir.getParent().getParent()

    def noRawImageSentinel = file("${repoRootDir}/assets/NO_RAW_IMAGE")
    def noProbZarrSentinel = file("${repoRootDir}/assets/NO_PROB_ZARR")
    def noTilesCsv         = file("${repoRootDir}/assets/NO_TILES_CSV")
    def noLabelsZarr       = file("${repoRootDir}/assets/NO_LABELS_ZARR")
    def noContoursZarr     = file("${repoRootDir}/assets/NO_CONTOURS_ZARR")

    ch_split_ultrack_input = ch_ultrack_input
        .map { meta, prob_zarr, raw_image ->
            def (localRawImage, rawImageUri) = Utils.splitImageInput(raw_image, noRawImageSentinel)
            def (localProbZarr, probZarrUri) = Utils.splitImageInput(prob_zarr, noProbZarrSentinel)
            tuple(meta.id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri)
        }

    ch_prealloc_input = ch_split_ultrack_input
        .map { id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri ->
            def prealloc_meta = meta.clone()
            prealloc_meta.num_channels = 1
            prealloc_meta.dtype        = "int32"
            prealloc_meta.data_type    = "label"
            prealloc_meta.zarr_output = meta.output_zarr
            tuple(prealloc_meta, localRawImage, rawImageUri, noTilesCsv)
        }

    PREALLOCATE_ZARR(ch_prealloc_input)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_ultrack = ch_split_ultrack_input
        .combine(
            PREALLOCATE_ZARR.out.zarr.map { meta, zarr -> tuple(meta.id, zarr) },
            by: 0
        )
        .map { id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri, prealloc_zarr ->
            tuple(
                meta,
                noLabelsZarr, null,
                localProbZarr, probZarrUri,
                noContoursZarr, null,
                localRawImage, rawImageUri,
                prealloc_zarr
            )
        }

    RUN_ULTRACK(ch_ultrack, val_config_toml)
    ch_versions = ch_versions.mix(RUN_ULTRACK.out.versions)

    ch_ready_for_pyramid = RUN_ULTRACK.out.tracked_zarr
        .map { meta, output_zarr -> tuple(meta, output_zarr, "Sample") }

    POPULATE_PYRAMID(ch_ready_for_pyramid)
    ch_versions = ch_versions.mix(POPULATE_PYRAMID.out.versions)

    emit:
    tracked_zarr = POPULATE_PYRAMID.out.zarr
    versions     = ch_versions
}
