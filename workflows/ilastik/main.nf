#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREALLOCATE_OUTPUT_ZARR  } from '../../modules/local/preallocate_output_zarr/main'
include { TILE_CONSTRUCTOR  } from '../../modules/local/tile_constructor/main'
include { ILASTIK_SUBREGION_PIXEL_CLASS  } from '../../modules/local/ilastik_segment/main'

workflow {
    // 1. Read parent images and metadata from the input CSV file
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

    // Define configuration values
    def ilastik_project = file(params.ilastik_project)
    def target_pubdir   = "${params.outdir}/probabilities"

    // 2. Preallocate the empty OME-Zarr skeleton on the shared disk
    PREALLOCATE_OUTPUT_ZARR(input_ch, target_pubdir) 

    // 3. Slice the parent images into coordinate tiles (outputs a tiles CSV)
    TILE_CONSTRUCTOR(input_ch)
    
    // 4. Parse the generated tiles CSV into individual task metadata
    ch_image_tiles = TILE_CONSTRUCTOR.out.tiles
        .map { _parentMeta, tileCsv -> tileCsv }
        .splitCsv(header: true)
        .map { row ->
            def tileMeta = [
                id: row.dataset_id,                       // e.g., "sample_01__tile_000001"
                source_dataset_id: row.source_dataset_id, // e.g., "sample_01" (joins with preallocated skeleton)
                resolution_level: row.resolution_level, 
                x_min: row.x_min, x_max: row.x_max, 
                y_min: row.y_min, y_max: row.y_max,
                z_min: row.z_min, z_max: row.z_max
            ]
            // We emit source_dataset_id ("sample_01") as index 0 to pair with the skeleton
            tuple(row.source_dataset_id, tileMeta, row.input_uri)
        }

    // 5. Prepare the preallocated skeleton channel for the join step
    ch_preallocated_ready = PREALLOCATE_OUTPUT_ZARR.out.empty_zarr
        .map { meta, zarr_path -> tuple(meta.id, zarr_path) }

    // 6. Join the individual tiles with their parent's preallocated Zarr store
    ch_joined_inputs = ch_image_tiles
        .join(ch_preallocated_ready, by: 0)
        .map { source_id, tileMeta, input_uri, preallocated_zarr_path ->
            // Format to match: [ val(tileMeta), val(input_uri), val(preallocated_zarr_path) ]
            tuple(tileMeta, input_uri, preallocated_zarr_path)
        }

    // 7. Run parallel pixel classification writing directly to the target Zarr
    ILASTIK_SUBREGION_PIXEL_CLASS(ch_joined_inputs, ilastik_project)
}
