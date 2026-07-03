#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SLICE_IMAGE  } from '../../modules/local/split_channels/main'
include { ILASTIK_SEGMENT  } from '../../modules/local/ilastik_segment/main'

workflow {
    // 1. Ingest raw CSV data rows
    input_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                sample: row.sample,
                exp_name: row.exp_name,
                csv_axes: row.axes,
                csv_source: row.export_source,
                csv_dtype: row.dtype
            ]
            return [ meta, file(row.file_path) ]
        }

    // 2. Split tif stack into separate frame slices
    slice_output = SLICE_IMAGE(input_ch)

    // 3. Flatten out the slice files array and build individual meta tracking components
    slice_ch = slice_output.out.single_slices
        .transpose()
        .map { meta, slice_file ->
            def matcher = (slice_file.name =~ /_t(\d+)_z(\d+)_c(\d+)\.tif$/)
            def slice_identity = matcher ? "t${matcher[0][1]}_z${matcher[0][2]}_c${matcher[0][3]}" : "unknown"

            return [
                meta + [
                    slice_id: slice_identity,
                    id: "${meta.sample}_${slice_identity}"
                ],
                slice_file
            ]
        }
    // 4. Run ilastik for each slice in parallel
    ILASTIK_SEGMENT(slice_ch, params.ilastik_project)
}
