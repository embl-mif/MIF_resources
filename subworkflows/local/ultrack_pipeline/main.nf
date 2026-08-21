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

    ch_prealloc_input = ch_ultrack_input
        .map { meta, prob_zarr, raw_image ->
            def prealloc_meta = meta.clone()
            prealloc_meta.num_channels = 1
            prealloc_meta.dtype        = "int32"
            prealloc_meta.data_type    = "label"
            def output_tracked_zarr = "${meta.output_dir}/${meta.id}_tracked.zarr"
            tuple(prealloc_meta, raw_image, output_tracked_zarr)
        }

    PREALLOCATE_ZARR(ch_prealloc_input)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_ultrack = ch_ultrack_input
        .map { meta, prob_zarr, raw_image -> tuple(meta.id, meta, prob_zarr, raw_image) }
        .combine(
            PREALLOCATE_ZARR.out.empty_zarr.map { meta, zarr -> tuple(meta.id, zarr) },
            by: 0
        )
        .map { id, meta, prob_zarr, raw_image, prealloc_zarr ->
            tuple(meta, null, prob_zarr, null, raw_image, prealloc_zarr)
        }

    RUN_ULTRACK(ch_ultrack, val_config_toml)
    ch_versions = ch_versions.mix(RUN_ULTRACK.out.versions)

    // downsample method = "Sample" for integer labels, not "Average"
    ch_ready_for_pyramid = RUN_ULTRACK.out.tracked_zarr
        .map { meta, output_zarr -> tuple(meta, output_zarr, "Sample") }

    POPULATE_PYRAMID(ch_ready_for_pyramid)
    ch_versions = ch_versions.mix(POPULATE_PYRAMID.out.versions)

    emit:
    tracked_zarr = POPULATE_PYRAMID.out.zarr
    versions     = ch_versions
}
