#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TILECONSTRUCTOR_INPUT_PREP    } from '../../modules/local/tileconstructor_input_prep/main'
include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { ILASTIK_ZARR }                  from '../../modules/local/ilastik_zarr/main'
include { REBUILD_PYRAMID               } from '../../modules/local/rebuild_pyramid/main'

workflow {
    def pubdir = "${file(params.outdir)}/probabilities"

    config_csv_ch = TILECONSTRUCTOR_INPUT_PREP(params.input_dir, params.resolution_level)

    // Read parent images and metadata from the input CSV file
    input_ch = config_csv_ch.csv
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id               : row.dataset_id,
                tile_size_x      : row.tile_size_x as int,
                tile_size_y      : row.tile_size_y as int,
                tile_size_z      : row.tile_size_z as int,
                tile_overlap     : row.tile_overlap as float,
                resolution_level : row.resolution_level as int,
                x_min            : row.x_min ? row.x_min as float : null,
                x_max            : row.x_max ? row.x_max as float : null,
                y_min            : row.y_min ? row.y_min as float : null,
                y_max            : row.y_max ? row.y_max as float : null,
                z_min            : row.z_min ? row.z_min as float : null,
                z_max            : row.z_max ? row.z_max as float : null,
                t_min            : row.t_min != null && row.t_min != '' ? row.t_min as int : null,
                t_max            : row.t_max != null && row.t_max != '' ? row.t_max as int : null,
                use_physical_units: false
            ]
            def image = file(row.input_image_path)
            return [ meta, image ]
        }

    // Wrap ilastik project file in a value channel so it can be reused by all parallel tasks
    ch_project = Channel.value(file(params.ilastik_project))

    TILE_CONSTRUCTOR(input_ch)

    ch_ilastik_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(input_ch)   // -> [ meta, tiles_csv, raw_zarr ]
        .flatMap { meta, tiles_csv, raw_zarr ->
            def output_zarr = "${pubdir}/${meta.id}_probabilities.zarr"

            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta.clone()
                tile_meta.y_min = row.y_min as Integer
                tile_meta.y_max = row.y_max as Integer
                tile_meta.x_min = row.x_min as Integer
                tile_meta.x_max = row.x_max as Integer
                tile_meta.z_min = (row.z_min && row.z_min != '') ? row.z_min as Integer : null
                tile_meta.z_max = (row.z_max && row.z_max != '') ? row.z_max as Integer : null
                // t_min/t_max unchanged per tile

                return [ tile_meta, raw_zarr, output_zarr ]
            }
        }

    ILASTIK_ZARR(ch_ilastik_inputs, ch_project)

    ch_ready_for_pyramid = ILASTIK_ZARR.out.zarr
        .map { meta, output_zarr -> tuple(meta.id, meta, output_zarr) }
        .groupTuple(by: 0)
        .map { id, metas, zarrs -> tuple(metas[0], zarrs[0], "Average") } // downsample method = average if using prob

    REBUILD_PYRAMID(ch_ready_for_pyramid)
}
