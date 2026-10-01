#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION as PYMIF_CONVERSION_INTENSITY         } from '../../modules/local/pymif_conversion/main'
include { PYMIF_CONVERSION as PYMIF_CONVERSION_LABELS } from '../../modules/local/pymif_conversion/main'
include { TILED_CELLPOSE_SAM_ZARR  } from '../../subworkflows/local/tiled_cellpose_sam_zarr/main'
include { PYMIF_FEATURE_EXTRACTION } from '../../modules/local/pymif_feature_extraction/main'

workflow {
    def default_features_dir = "${file(params.outdir)}/measurements"

    // Pipeline-wide, not per-image -- label_offset_step must be large enough
    // that no two tiles' own label IDs can collide once offset (tile_index *
    // label_offset_step), and multicut_beta is LABEL_SOLVER's own edge-cut
    // trade-off parameter (label_solver.py itself defaults to 0.5 if omitted).
    def label_offset_step = params.label_offset_step ?: 100000
    def multicut_beta     = params.multicut_beta     ?: 0.5

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()

            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            // cellpose tiling parameters
            def tile_size_x  = Utils.normalizeOptionalText(row.tile_size_x)
            def tile_size_y  = Utils.normalizeOptionalText(row.tile_size_y)
            def tile_size_z  = Utils.normalizeOptionalText(row.tile_size_z)
            def tile_overlap = Utils.normalizeOptionalText(row.tile_overlap)

            def resolution_level   = Utils.normalizeOptionalText(row.resolution_level)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)

            def x_min = Utils.normalizeOptionalText(row.x_min)
            def x_max = Utils.normalizeOptionalText(row.x_max)
            def y_min = Utils.normalizeOptionalText(row.y_min)
            def y_max = Utils.normalizeOptionalText(row.y_max)
            def z_min = Utils.normalizeOptionalText(row.z_min)
            def z_max = Utils.normalizeOptionalText(row.z_max)

            def diameter_value    = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value       = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo     = Utils.normalizeOptionalText(row.cellpose_halo)
            def flow_threshold    = Utils.normalizeOptionalText(row.cellpose_flow_threshold)
            def cellprob_threshold = Utils.normalizeOptionalText(row.cellpose_cellprob_threshold)

            def do_3D_text  = Utils.normalizeOptionalText(row.do_3D)
            def do_3D_value = do_3D_text ? do_3D_text.equalsIgnoreCase('true') : false

            def resolved_zarr_output = Utils.resolveOutput(row, 'output', "${dataset_id}.ome.zarr", dataset_id)
            def resolved_labels_output = Utils.resolveOutput(row, 'labels_output_final', "${dataset_id}_merged_mask.ome.zarr", dataset_id)

            // Feature extraction parameters
            def intensity_resolution_level = Utils.normalizeOptionalText(row.intensity_resolution_level)
            def labels_resolution_level = Utils.normalizeOptionalText(row.labels_resolution_level)

            def resolved_features_output = Utils.resolveOutput(row, 'features_output', "${default_features_dir}/${dataset_id}_features.csv", dataset_id)

            def meta = row + [
                id                 : dataset_id,
                zarr_output_name   : resolved_zarr_output,
                labels_output_final: resolved_labels_output,
                // The subworkflow writes its merged mask to meta.labels_output. Point it at the
                // work dir so it stays intermediate; the CSV value is kept in labels_output_final.
                labels_output      : "${dataset_id}_merged_mask.ome.zarr",
                features_output_name: resolved_features_output,

                tile_size_x        : tile_size_x  ? tile_size_x  as int   : null,
                tile_size_y        : tile_size_y  ? tile_size_y  as int   : null,
                tile_size_z        : tile_size_z  ? tile_size_z  as int   : null,
                tile_overlap       : tile_overlap ? tile_overlap as float : null,
                resolution_level   : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : null),
                use_physical_units : use_physical_units ? use_physical_units.toBoolean() : false,

                x_min : x_min ? x_min as float : null,
                x_max : x_max ? x_max as float : null,
                y_min : y_min ? y_min as float : null,
                y_max : y_max ? y_max as float : null,
                z_min : z_min ? z_min as float : null,
                z_max : z_max ? z_max as float : null,

                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                flow_threshold        : flow_threshold ? flow_threshold as float : null,
                cellprob_threshold    : cellprob_threshold ? cellprob_threshold as float : null,
                halo : cellpose_halo != null && cellpose_halo != '' ? cellpose_halo as float : 20.0,
                do_3D : do_3D_value,

                features_intensity_resolution_level       : intensity_resolution_level
                    ? intensity_resolution_level as int
                    : (params.intensity_resolution_level != null ? params.intensity_resolution_level as int : 0),
                features_labels_resolution_level          : labels_resolution_level
                    ? labels_resolution_level as int
                    : (params.labels_resolution_level != null ? params.labels_resolution_level as int : 0),
                features_channels                         : Utils.normalizeOptionalText(row.channels),
                features_extra_properties                 : Utils.normalizeOptionalText(row.extra_properties),
                features_extra_intensity_properties       : Utils.normalizeOptionalText(row.extra_intensity_properties),

            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION_INTENSITY(pymif_inputs_ch)

    ch_images = PYMIF_CONVERSION_INTENSITY.out.zarr

    TILED_CELLPOSE_SAM_ZARR(
        ch_images,
        label_offset_step,
        multicut_beta,
    )

    TILED_CELLPOSE_SAM_ZARR.out.mask.view { meta, mask_path ->
        "[TILED_CELLPOSE_SAM_ZARR] ${meta.id}: intermediate merged mask ready at ${mask_path}"
    }

    // Convert the single-resolution intermediate mask into the final label OME-Zarr with
    // pymif (zarr -> zarr), reusing the row's conversion settings (num_levels,
    // downscale_factor, zarr_format, ...) so the pyramid matches the intensity image.
    // Intensity-only options (channels, chunk_size, subset) don't apply to the mask.
    labels_inputs_ch = TILED_CELLPOSE_SAM_ZARR.out.mask
        .map { meta, mask -> tuple(meta.id, mask) }
        .join(ch_images.map { meta, _zarr -> tuple(meta.id, meta) })
        .map { _id, mask, meta ->
            def labels_meta = meta + [
                microscope      : 'zarr',
                zarr_output_name: meta.labels_output_final,
                chunk_size      : null,
                channel_names   : null,
                channel_colors  : null,
                subset          : null,
                scene_index     : null,
            ]
            tuple(labels_meta, file(mask.toString()))
        }

    PYMIF_CONVERSION_LABELS(labels_inputs_ch)

    PYMIF_CONVERSION_LABELS.out.zarr.view { meta, labels_path ->
        "[PYMIF_CONVERSION_LABELS] ${meta.id}: label OME-Zarr ready at ${labels_path}"
    }

    // Measure the converted image against the final label OME-Zarr.
    feature_inputs_ch = ch_images
        .map { meta, zarr -> tuple(meta.id, meta, zarr) }
        .join(PYMIF_CONVERSION_LABELS.out.zarr.map { meta, labels -> tuple(meta.id, labels) })
        .map { _id, meta, intensity, labels -> tuple(meta, intensity, labels) }

    PYMIF_FEATURE_EXTRACTION(feature_inputs_ch)

}
