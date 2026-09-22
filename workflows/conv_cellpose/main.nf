#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'


workflow {

    def resolveOutput = { value, fallback ->
        def v = Utils.normalizeOptionalText(value)
        !v ? fallback : (v.startsWith('/') ? v : file(v).name)
    }

    def repoRootDir = workflow.projectDir.getParent().getParent()

    def default_zarr_dir = "${file(params.default_outdir)}/zarr"
    def default_labels_dir = "${file(params.default_outdir)}/labels"

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${repoRootDir}/assets/NO_MODEL_FILE")
    )

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

            def diameter_value = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo = Utils.normalizeOptionalText(row.cellpose_halo)
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)
            def flow_threshold     = Utils.normalizeOptionalText(row.cellpose_flow_threshold)
            def cellprob_threshold = Utils.normalizeOptionalText(row.cellpose_cellprob_threshold)

            def x_min = Utils.normalizeOptionalText(row.x_min)
            def x_max = Utils.normalizeOptionalText(row.x_max)
            def y_min = Utils.normalizeOptionalText(row.y_min)
            def y_max = Utils.normalizeOptionalText(row.y_max)
            def z_min = Utils.normalizeOptionalText(row.z_min)
            def z_max = Utils.normalizeOptionalText(row.z_max)

            def do_3D_text = Utils.normalizeOptionalText(row.do_3D)
            def do_3D_value = do_3D_text
                ? do_3D_text.equalsIgnoreCase('true')
                : false

            def requested_zarr_output = Utils.normalizeOptionalText(row.output)
            if (requested_zarr_output && !file(requested_zarr_output).isAbsolute()) {
                error(
                    "zarr_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_zarr_output}'. Leave the column blank to use the " +
                    "default location under outdir instead."
                )
            }
            def resolved_zarr_output = requested_zarr_output
                ?: "${default_zarr_dir}/${dataset_id}.ome.zarr"


            def requested_labels_output = Utils.normalizeOptionalText(row.labels_output)
            if (requested_labels_output && !file(requested_labels_output).isAbsolute()) {
                error(
                    "labels_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_labels_output}'. Leave the column blank to use the " +
                    "default location under outdir instead."
                )
            }
            def resolved_labels_output = requested_labels_output
                ?: "${default_labels_dir}/${dataset_id}_labels.ome.zarr"

            def meta = row + [
                id                    : dataset_id,
                zarr_output_name      : resolved_zarr_output,
                labels_output_name    : resolved_labels_output,
                resolution_level      : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : null),
                use_physical_units    : use_physical_units ? use_physical_units.toBoolean() : false,
                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                flow_threshold        : flow_threshold ? flow_threshold as float : null,
                cellprob_threshold    : cellprob_threshold ? cellprob_threshold as float : null,
                halo                  : cellpose_halo != null && cellpose_halo != ''
                    ? cellpose_halo as float
                    : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                x_min : x_min ? x_min as float : null,
                x_max : x_max ? x_max as float : null,
                y_min : y_min ? y_min as float : null,
                y_max : y_max ? y_max as float : null,
                z_min : z_min ? z_min as float : null,
                z_max : z_max ? z_max as float : null,
                do_3D                 : do_3D_value
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    ch_cellpose_inputs = PYMIF_CONVERSION.out.zarr
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${repoRootDir}/assets/NO_FILE")
            )
            tuple(meta, local_image, image_uri)
        }
 
    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)
}

