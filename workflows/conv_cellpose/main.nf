#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

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

            def diameter_value = row.cellpose_diameter?.toString()?.trim()
            def niter_value = row.cellpose_niter?.toString()?.trim()
            def cellpose_channels = row.cellpose_channels?.toString()?.trim()

            def resolved_zarr_output =
                "${zarr_dir}/${file(row.output).name}"
            
            def labels_output_value = row.labels_output?.toString()?.trim()
            def labels_output_name = labels_output_value
                ? file(labels_output_value).name
                : "${dataset_id}_labels.ome.zarr"
            def resolved_labels_output = "${zarr_dir}/${labels_output_name}"

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,
                labels_output_name       : resolved_labels_output,
                diameter                 : diameter_value ? diameter_value as int : null,
                niter                    : niter_value ? niter_value as int : null,
                segmentation_channels    : cellpose_channels ?: null
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    CELLPOSE_SAM_ZARR_MIF(PYMIF_CONVERSION.out.zarr)
}

