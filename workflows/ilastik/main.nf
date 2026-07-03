#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SLICE_IMAGE  } from '../../modules/local/split_channels/main'
include { ILASTIK_SEGMENT  } from '../../modules/local/ilastik_segment/main'
include { HYPERSTACK  } from '../../modules/local/hyperstack/main'
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
                csv_dtype: row.dtype,
                target_channel: row.target_channel ? row.target_channel.toInteger() : 0
            ]
            return [ meta, file(row.file_path) ]
        }

    // 2. Split TCZYX OME-TIFF hyperstack into separate frame slices
    slice_output = SLICE_IMAGE(input_ch)

    // 3. Flatten out the slice files array and build individual meta tracking components
    slice_ch = slice_output.out.single_slices
        .transpose()
        .map { meta, slice_file ->
            def matcher = (slice_file.name =~ /_t(\d+)_z(\d+)_c(\d+)\.tif$/)
            def channel_string = matcher ? matcher[0][3] : "${meta.target_channel}"
            def slice_identity = matcher ? "t${matcher[0][1]}_z${matcher[0][2]}_c${channel_string}" : "unknown"

            return [
                meta + [
                    slice_id: slice_identity,
                    channel_idx: channel_string.toInteger(),
                    id: "${meta.sample}_${slice_identity}"
                ],
                slice_file
            ]
        }
        .filter { meta, slice_file -> 
            meta.channel_idx == meta.target_channel
        }
    // 4. Run ilastik for each slice in parallel
    ilastik_out = ILASTIK_SEGMENT(slice_ch, params.ilastik_project)

    // 5. Re-group individual slices by sample name
    grouped_masks_ch = ilastik_out.single_segs
        .map { meta, seg_file ->
            // Revert 'id' back to original image sample name
            def parent_meta = [
                id       : meta.sample,
                sample   : meta.sample,
                exp_name : meta.exp_name
            ]
            return [ parent_meta, seg_file ]
        }
        .groupTuple(by: 0)

    // 6. Hyperstack
    HYPERSTACK(grouped_masks_ch)
}
