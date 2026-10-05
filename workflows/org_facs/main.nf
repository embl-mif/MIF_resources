#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION         } from '../../modules/local/pymif_conversion/main'
include { TILED_CELLPOSE_SAM_ZARR  } from '../../subworkflows/local/tiled_cellpose_sam_zarr/main'
include { PYMIF_FEATURE_EXTRACTION } from '../../modules/local/pymif_feature_extraction/main'

workflow {
    def default_features_dir = "${file(params.outdir)}/measurements"
    def default_labels_dir   = "${file(params.outdir)}/masks"

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

            def zarr_output = Utils.resolveOutput(row, 'output', "${dataset_id}.ome.zarr", dataset_id)

            tuple(row + [id: dataset_id, zarr_output: zarr_output], input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    ch_images = PYMIF_CONVERSION.out.zarr

    def withCellposeSettings = { meta ->
        def tile_size_x  = Utils.normalizeOptionalText(meta.tile_size_x)
        def tile_size_y  = Utils.normalizeOptionalText(meta.tile_size_y)
        def tile_size_z  = Utils.normalizeOptionalText(meta.tile_size_z)
        def tile_overlap = Utils.normalizeOptionalText(meta.tile_overlap)

        def resolution_level   = Utils.normalizeOptionalText(meta.resolution_level)
        def use_physical_units = Utils.normalizeOptionalText(meta.use_physical_units)

        def x_min = Utils.normalizeOptionalText(meta.x_min)
        def x_max = Utils.normalizeOptionalText(meta.x_max)
        def y_min = Utils.normalizeOptionalText(meta.y_min)
        def y_max = Utils.normalizeOptionalText(meta.y_max)
        def z_min = Utils.normalizeOptionalText(meta.z_min)
        def z_max = Utils.normalizeOptionalText(meta.z_max)

        def diameter_value     = Utils.normalizeOptionalText(meta.cellpose_diameter)
        def niter_value        = Utils.normalizeOptionalText(meta.cellpose_niter)
        def cellpose_channels  = Utils.normalizeOptionalText(meta.cellpose_channels)
        def cellpose_halo      = Utils.normalizeOptionalText(meta.cellpose_halo)
        def flow_threshold     = Utils.normalizeOptionalText(meta.cellpose_flow_threshold)
        def cellprob_threshold = Utils.normalizeOptionalText(meta.cellpose_cellprob_threshold)

        def do_3D_text  = Utils.normalizeOptionalText(meta.do_3D)
        def do_3D_value = do_3D_text ? do_3D_text.equalsIgnoreCase('true') : false

        meta + [
            labels_output      : Utils.resolveOutput(meta, 'labels_output_final', "${default_labels_dir}/${meta.id}_merged_mask.ome.zarr", meta.id),

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
        ]
    }

    TILED_CELLPOSE_SAM_ZARR(
        ch_images.map { meta, zarr -> tuple(withCellposeSettings(meta), zarr) },
        label_offset_step,
        multicut_beta,
    )

    TILED_CELLPOSE_SAM_ZARR.out.labels.view { meta, mask_path ->
        "[TILED_CELLPOSE_SAM_ZARR] ${meta.id}: final label OME-Zarr ready at ${mask_path}"
    }

    feature_inputs_ch = ch_images
        .map { meta, zarr -> tuple(meta.id, meta, zarr) }
        .join(TILED_CELLPOSE_SAM_ZARR.out.labels.map { meta, labels -> tuple(meta.id, labels) })
        .map { _id, meta, intensity, labels ->
            def intensity_resolution_level = Utils.normalizeOptionalText(meta.intensity_resolution_level)
            def labels_resolution_level    = Utils.normalizeOptionalText(meta.labels_resolution_level)

            def features_meta = meta + [
                csv_output: Utils.resolveOutput(meta, 'features_output', "${default_features_dir}/${meta.id}_features.csv", meta.id),

                features_intensity_resolution_level       : intensity_resolution_level
                    ? intensity_resolution_level as int
                    : (params.intensity_resolution_level != null ? params.intensity_resolution_level as int : 0),
                features_labels_resolution_level          : labels_resolution_level
                    ? labels_resolution_level as int
                    : (params.labels_resolution_level != null ? params.labels_resolution_level as int : 0),
                features_channels                         : Utils.normalizeOptionalText(meta.channels),
                features_extra_properties                 : Utils.normalizeOptionalText(meta.extra_properties),
                features_extra_intensity_properties       : Utils.normalizeOptionalText(meta.extra_intensity_properties),
            ]
            tuple(features_meta, intensity, labels)
        }

    PYMIF_FEATURE_EXTRACTION(feature_inputs_ch)

}
