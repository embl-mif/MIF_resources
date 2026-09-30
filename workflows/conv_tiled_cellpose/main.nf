#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION         } from '../../modules/local/pymif_conversion/main'
include { TILED_CELLPOSE_SAM_ZARR  } from '../../subworkflows/local/tiled_cellpose_sam_zarr/main'

workflow {
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

            // A relative `input` value is resolved against the CSV's own
            // directory (not wherever Nextflow was launched from) -- see
            // Utils.resolvePath. Already returns an absolute, existence-
            // checked Path for the local-file case.
            def input_dataset = Utils.resolvePath(row.input, params.input_csv)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

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

            // Whole-workflow policy: an output column left blank in the CSV
            // defaults to a bare relative filename, so it simply lives in
            // (and is cleaned up along with) the task's own work directory --
            // nothing is written to persistent/shared storage unless the CSV
            // explicitly gives an absolute path for that column.
            def resolved_zarr_output = Utils.resolveOutput(row, 'output', "${dataset_id}.ome.zarr", dataset_id)
            def resolved_labels_output = Utils.resolveOutput(row, 'labels_output', "${dataset_id}_merged_mask.ome.zarr", dataset_id)

            def meta = row + [
                id                 : dataset_id,
                zarr_output_name   : resolved_zarr_output,
                labels_output      : resolved_labels_output,

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

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    ch_images = PYMIF_CONVERSION.out.zarr

    TILED_CELLPOSE_SAM_ZARR(
        ch_images,
        label_offset_step,
        multicut_beta,
    )

    TILED_CELLPOSE_SAM_ZARR.out.mask.view { meta, mask_path ->
        "[TILED_CELLPOSE_SAM_ZARR] ${meta.id}: merged mask ready at ${mask_path}"
    }
}
