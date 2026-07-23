#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TILECONSTRUCTOR_INPUT_PREP  } from '../../modules/local/tileconstructor_input_prep/main'
include { PREALLOCATE_OUTPUT_ZARR  } from '../../modules/local/preallocate_output_zarr/main'
include { TILE_CONSTRUCTOR  } from '../../modules/local/tileconstructor/main'
include { ILASTIK_SUBREGION_PIXEL_CLASS  } from '../../modules/local/ilastik_subregion_pixel_class/main'

workflow {
    def pubdir = "${file(params.outdir)}/probabilities"

    config_csv_ch = TILECONSTRUCTOR_INPUT_PREP(params.input_dir, params.resolution_level)

    //  Read parent images and metadata from the input CSV file
    input_ch = config_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id            : row.dataset_id,
                tile_size_x   : row.tile_size_x as int,
                tile_size_y   : row.tile_size_y as int,
                tile_size_z   : row.tile_size_z as int,
                tile_overlap  : row.tile_overlap as float,
                resolution_level: row.resolution_level as int,
                x_min         : row.x_min ? row.x_min as float : null,
                x_max         : row.x_max ? row.x_max as float : null,
                y_min         : row.y_min ? row.y_min as float : null,
                y_max         : row.y_max ? row.y_max as float : null,
                z_min         : row.z_min ? row.z_min as float : null,
                z_max         : row.z_max ? row.z_max as float : null,
                use_physical_units: false
            ]
            def image = file(row.input_image_path)
            return [ meta, image ]
        }

    // Wrap ilastik project file in a value channel so it can be reused by all parallel tasks
    ch_project = Channel.value(file(params.ilastik_project))

    // Preallocate the empty OME-Zarr skeleton
    PREALLOCATE_OUTPUT_ZARR(input_ch, ch_project, pubdir)

    // Slice the parent images into coordinate tiles (outputs a tiles CSV)
    TILE_CONSTRUCTOR(input_ch)
    
    ch_ilastik_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(input_ch)                      // Combines -> [ meta, tiles_csv, raw_zarr ]
        .join(PREALLOCATE_OUTPUT_ZARR.out.empty_zarr)   // Combines -> [ meta, tiles_csv, raw_zarr, preallocated_zarr ]
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
                tile_meta.z_index = (row.z_index != null && row.z_index != '') ? row.z_index as Integer : null
                tile_meta.t_index = (row.t_index != null && row.t_index != '') ? row.t_index as Integer : (row.t != null && row.t != '' ? row.t as Integer : 0)
                return [ tile_meta, raw_zarr, preallocated_zarr ]
            }
        }
    ILASTIK_SUBREGION_PIXEL_CLASS(ch_ilastik_inputs, ch_project)
}
