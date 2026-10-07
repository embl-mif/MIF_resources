#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { ILASTIK_ZARR_PIPELINE         } from '../../subworkflows/local/ilastik_zarr_pipeline/main'

workflow {
    def pixel_class_dir = "${file(params.outdir)}/pixel_class"

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    input_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()

            def input_dataset = file(row.input_image_path, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def tile_size_x  = Utils.normalizeOptionalText(row.tile_size_x)
            def tile_size_y  = Utils.normalizeOptionalText(row.tile_size_y)
            def tile_size_z  = Utils.normalizeOptionalText(row.tile_size_z)
            def tile_overlap = Utils.normalizeOptionalText(row.tile_overlap)
            def halo         = Utils.normalizeOptionalText(row.halo)

            def resolution_level   = Utils.normalizeOptionalText(row.resolution_level)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)

            def x_min = Utils.normalizeOptionalText(row.x_min)
            def x_max = Utils.normalizeOptionalText(row.x_max)
            def y_min = Utils.normalizeOptionalText(row.y_min)
            def y_max = Utils.normalizeOptionalText(row.y_max)
            def z_min = Utils.normalizeOptionalText(row.z_min)
            def z_max = Utils.normalizeOptionalText(row.z_max)

            def t_min = Utils.normalizeOptionalText(row.t_min)
            def t_max = Utils.normalizeOptionalText(row.t_max)

            def resolved_output_zarr = Utils.resolveOutput(row, 'pixel_class_output', "${pixel_class_dir}/${dataset_id}_pixel_class.zarr", dataset_id)

            def meta = row + [
                id                 : dataset_id,
                zarr_output        : resolved_output_zarr,

                tile_size_x        : tile_size_x  ? tile_size_x  as int   : null,
                tile_size_y        : tile_size_y  ? tile_size_y  as int   : null,
                tile_size_z        : tile_size_z  ? tile_size_z  as int   : null,
                tile_overlap       : tile_overlap ? tile_overlap as float : null,
                halo               : halo ? halo as float : null,
                resolution_level   : resolution_level ? resolution_level as int : null,
                use_physical_units : use_physical_units ? use_physical_units.toBoolean() : false,

                x_min : x_min ? x_min as float : null,
                x_max : x_max ? x_max as float : null,
                y_min : y_min ? y_min as float : null,
                y_max : y_max ? y_max as float : null,
                z_min : z_min ? z_min as float : null,
                z_max : z_max ? z_max as float : null,

                t_min : t_min ? t_min as int : null,
                t_max : t_max ? t_max as int : null,
            ]

            tuple(meta, input_dataset)
        }

    TILE_CONSTRUCTOR(input_ch)

    ch_ilastik_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(input_ch)
        .flatMap { meta, tiles_csv, raw_zarr ->
            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta.clone()
                tile_meta.resolution_level = row.resolution_level as Integer
                tile_meta.y_min = row.y_min as Integer
                tile_meta.y_max = row.y_max as Integer
                tile_meta.x_min = row.x_min as Integer
                tile_meta.x_max = row.x_max as Integer
                tile_meta.z_min = (row.z_min && row.z_min != '') ? row.z_min as Integer : null
                tile_meta.z_max = (row.z_max && row.z_max != '') ? row.z_max as Integer : null
                // t_min/t_max unchanged per tile

                return [ tile_meta, raw_zarr, meta.zarr_output ]
            }
        }
    ch_project = Channel.value(file(params.ilastik_project))

    ILASTIK_ZARR_PIPELINE(ch_ilastik_inputs, ch_project)
    ILASTIK_ZARR_PIPELINE.out.probabilities.view()
}
