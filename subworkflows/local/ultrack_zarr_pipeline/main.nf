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

    def intOr = { meta, column, fallback ->
        def value = Utils.parseOptionalInt(meta[column], column)
        value != null ? value : fallback
    }
    def floatOr = { meta, column, fallback ->
        def value = Utils.parseOptionalFloat(meta[column], column)
        value != null ? value : fallback
    }

    def preallocMeta = { meta ->
        [
            id           : meta.id,
            num_channels : 1,
            dtype        : 'int32',
            data_type    : 'label',
            zarr_output  : meta.zarr_output,
            t_min        : intOr(meta, 't_min', 0),
            t_max        : intOr(meta, 't_max', 5),
        ]
    }

    def ultrackMeta = { meta ->
        [
            id                 : meta.id,
            zarr_output        : meta.zarr_output,
            foreground_channel : intOr(meta, 'foreground_channel', 1),
            contours_channel   : intOr(meta, 'contours_channel', null),
            raw_channel        : intOr(meta, 'raw_channel', 0),
            t_min              : intOr(meta, 't_min', 0),
            t_max              : intOr(meta, 't_max', 5),
            min_area           : intOr(meta, 'min_area', 200),
            max_area           : intOr(meta, 'max_area', 20000),
            max_distance       : floatOr(meta, 'max_distance', 10.0),
            distance_weight    : floatOr(meta, 'distance_weight', null),
            solution_gap       : floatOr(meta, 'solution_gap', 0.1),
            time_limit         : intOr(meta, 'time_limit', 600),
            scale_z            : floatOr(meta, 'scale_z', null),
            scale_y            : floatOr(meta, 'scale_y', null),
            scale_x            : floatOr(meta, 'scale_x', null),
        ]
    }

    ch_split_ultrack_input = ch_ultrack_input
        .map { meta, prob_zarr, raw_image ->
            def (localRawImage, rawImageUri) = Utils.splitImageInput(raw_image, noRawImageSentinel)
            def (localProbZarr, probZarrUri) = Utils.splitImageInput(prob_zarr, noProbZarrSentinel)
            tuple(meta.id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri)
        }

    ch_prealloc_input = ch_split_ultrack_input
        .map { id, meta, localRawImage, rawImageUri, localProbZarr, probZarrUri ->
            tuple(preallocMeta(meta), localRawImage, rawImageUri, noTilesCsv)
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
                ultrackMeta(meta),
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
