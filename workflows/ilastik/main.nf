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
                sample: row.sample,
                exp_name: row.exp_name,
                csv_axes: row.axes,
                csv_source: row.export_source,
                csv_dtype: row.dtype
            ]
            return [ meta, file(row.file_path) ]
        }

    // 2. Split tif stack into separate frame slices
    SLICE_IMAGE(input_ch)

    // 3. Flatten out the slice files array and build individual meta tracking components
    slice_ch = SLICE_IMAGE(input_ch).out.single_slices
        .transpose()
        .map { meta, slice_file ->
            // Pull the raw filename identifier (e.g. 't001_c001')
            def slice_identity = slice_file.baseName.replaceAll("sample_${meta.sample}_", "")

            // Merge original fields into schema mapping meta.id
            def updated_meta = meta + [ 
                slice_id: slice_identity,
                id: "${meta.sample}_${slice_identity}"
            ]

            return [ updated_meta, slice_file ]
        }

    // 4. Run ilastik for each slice in parallel
    ILASTIK_SEGMENT(slice_ch, params.ilastik_project)
}
