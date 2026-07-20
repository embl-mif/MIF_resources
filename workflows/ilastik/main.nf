#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREALLOCATE_OUTPUT_ZARR  } from '../../modules/local/preallocate_output_zarr/main'
include { TILE_CONSTRUCTOR  } from '../../modules/local/tileconstructor/main'
include { ILASTIK_SUBREGION_PIXEL_CLASS  } from '../../modules/local/ilastik_subregion_pixel_class/main'

workflow {
    def pubdir = "${file(params.outdir)}/probabilities"

    //  Read parent images and metadata from the input CSV file
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
                tile_overlap: 0,
            ]
            def input_zarr = file(row.file_path, checkIfExists: true).toString()
            return [ meta, input_zarr ]
        }

    // Wrap ilastik project file in a value channel so it can be reused by all parallel tasks
    ch_project = Channel.value(file(params.ilastik_project))

    // Preallocate the empty OME-Zarr skeleton
    PREALLOCATE_OUTPUT_ZARR(input_ch, pubdir) 

    // Slice the parent images into coordinate tiles (outputs a tiles CSV)
    TILE_CONSTRUCTOR(input_ch)
    
    ch_ilastik_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(input_ch)                      // Combines -> [ meta, tiles_csv, raw_zarr ]
        .join(PREALLOCATE_OUTPUT_ZARR.out)   // Combines -> [ meta, tiles_csv, raw_zarr, preallocated_zarr ]
        .flatMap { meta, tiles_csv, raw_zarr, preallocated_zarr ->
            // Parse CSV rows inside flatMap to emit individual tile tasks
            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta.clone()
                tile_meta.y_min = row.y_min as Integer
                tile_meta.y_max = row.y_max as Integer
                tile_meta.x_min = row.x_min as Integer
                tile_meta.x_max = row.x_max as Integer

                // Handle optional Z bounds for 2D datasets
                tile_meta.z_min = (row.z_min && row.z_min != '') ? row.z_min as Integer : null
                tile_meta.z_max = (row.z_max && row.z_max != '') ? row.z_max as Integer : null

                return [ tile_meta, raw_zarr, preallocated_zarr ]
            }
        }
    ILASTIK_SUBREGION_PIXEL_CLASS(ch_ilastik_inputs, ch_project)
}
