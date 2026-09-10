#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP_ZARR                } from '../../modules/local/pymif_mip_zarr/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {
    def outdir = "${file(params.outdir)}"
    // if (!params.cellpose_direct_output_dir) {
    //     params.cellpose_direct_output_dir = "${file(params.outdir)}/labels"
    // }

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${workflow.projectDir}/assets/NO_MODEL_FILE")
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

            def resolved_zarr_output =
                "${outdir}/${file(row.output).name}"

            def labels_mip_value = Utils.normalizeOptionalText(row.mip_output)
            def labels_mip_name = labels_mip_value
                ? file(labels_mip_value).name
                : "${dataset_id}_mip.ome.zarr"
            def resolved_mip_output = "${outdir}/${labels_mip_name}"

            def diameter_value = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo = Utils.normalizeOptionalText(row.cellpose_halo)
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)
            def do_3D_text = Utils.normalizeOptionalText(row.do_3D)
            def do_3D_value = do_3D_text
                ? do_3D_text.equalsIgnoreCase('true')
                : false
            def labels_output_value = Utils.normalizeOptionalText(row.labels_output)
            def resolved_labels_output = labels_output_value
                ? (file(labels_output_value).isAbsolute()
                    ? labels_output_value                    // abs path: keep as given
                    : file(labels_output_value).name)        // relative: strip to bare filename
                : "${dataset_id}_labels.ome.zarr"             // not provided: default

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                mip_output_name          : resolved_mip_output,

                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                halo                  : cellpose_halo != null && cellpose_halo != ''
                    ? cellpose_halo as float
                    : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                resolution_level      : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : null),
                use_physical_units    : use_physical_units ? use_physical_units.toBoolean() : false,
                do_3D                 : do_3D_value,
                labels_output_name    : resolved_labels_output
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    /*
     * Creates a multiscale 2D OME-Zarr by applying max() over Z at every
     * resolution level.
     *
     * Expected output:
     * tuple(meta, mip_zarr)
     */
    PYMIF_MIP_ZARR(PYMIF_CONVERSION.out.zarr)

    ch_cellpose_inputs = PYMIF_MIP_ZARR.out.zarr
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${workflow.projectDir}/assets/NO_FILE")
            )
            tuple(meta, local_image, image_uri)
        }

    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)
}

