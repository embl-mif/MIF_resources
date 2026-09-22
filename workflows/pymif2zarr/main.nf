#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PYMIF_CONVERSION } from '../../modules/local/pymif_conversion/main'


workflow {
    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    def default_zarr_dir = "${file(params.outdir)}/zarr"

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"
            
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
            
            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)
}
