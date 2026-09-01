#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION_WIP              } from '../../modules/local/pymif_conversion_wip/main'
include { PYMIF_MIP_ZARR_WIP                } from '../../modules/local/pymif_mip_zarr_wip/main'
include { CELLPOSE_SAM_ZARR_WIP         } from '../../modules/local/cellpose_sam_zarr_wip/main'

workflow {
    def zarr_dir = "${file(params.outdir)}/zarr"

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = row.dataset_id?.toString()?.trim()
                ?: input_dataset.simpleName

            def mip_output_value = row.mip_output?.toString()?.trim()
                ?: "${dataset_id}_mip.ome.zarr"

            def diameter_value = row.cellpose_diameter?.toString()?.trim()
            def niter_value = row.cellpose_niter?.toString()?.trim()
            def cellpose_channels = row.cellpose_channels?.toString()?.trim()

            def resolved_zarr_output =
                "${zarr_dir}/${file(row.output).name}"
            
            def resolved_mip_output =
                "${zarr_dir}/${mip_output_value}"

            def meta = row + [
                id      : dataset_id,
                zarr_output_name : resolved_zarr_output,
                mip_output_name  : resolved_mip_output,
                diameter         : diameter_value ? diameter_value as int : null,
                niter            : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION_WIP(pymif_inputs_ch)

    /*
     * Creates a multiscale 2D OME-Zarr by applying max() over Z at every
     * resolution level.
     *
     * Expected output:
     * tuple(meta, mip_zarr)
     */
    PYMIF_MIP_ZARR_WIP(PYMIF_CONVERSION_WIP.out.zarr)

    CELLPOSE_SAM_ZARR_WIP(PYMIF_MIP_ZARR_WIP.out.zarr)

}

