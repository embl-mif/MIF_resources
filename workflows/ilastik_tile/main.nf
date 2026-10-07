#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { ILASTIK_ZARR_PIPELINE         } from '../../subworkflows/local/ilastik_zarr_pipeline/main'

workflow {
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

            def default_zarr_output = params.outdir
                ? "${file(params.outdir)}/pixel_class/${dataset_id}_pixel_class.zarr"
                : null

            if (!Utils.normalizeOptionalText(row.pixel_class_output) && !default_zarr_output) {
                throw new IllegalArgumentException(
                    "pixel_class_output is blank for dataset_id=${dataset_id} and --outdir is not set. " +
                    "Set an absolute pixel_class_output in the CSV or pass --outdir."
                )
            }

            def zarr_output = Utils.resolveOutput(row, 'pixel_class_output', default_zarr_output, dataset_id)

            tuple(row + [id: dataset_id, zarr_output: zarr_output], input_dataset)
        }

    def tileConstructorMeta = { meta ->
        def tile_size_x  = Utils.normalizeOptionalText(meta.tile_size_x)
        def tile_size_y  = Utils.normalizeOptionalText(meta.tile_size_y)
        def tile_size_z  = Utils.normalizeOptionalText(meta.tile_size_z)
        def tile_overlap = Utils.normalizeOptionalText(meta.tile_overlap)
        def halo         = Utils.normalizeOptionalText(meta.halo)

        def resolution_level   = Utils.normalizeOptionalText(meta.resolution_level)
        def use_physical_units = Utils.normalizeOptionalText(meta.use_physical_units)

        def x_min = Utils.normalizeOptionalText(meta.x_min)
        def x_max = Utils.normalizeOptionalText(meta.x_max)
        def y_min = Utils.normalizeOptionalText(meta.y_min)
        def y_max = Utils.normalizeOptionalText(meta.y_max)
        def z_min = Utils.normalizeOptionalText(meta.z_min)
        def z_max = Utils.normalizeOptionalText(meta.z_max)

        [
            id                 : meta.id,

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
        ]
    }

    def ilastikMeta = { meta ->
        def halo  = Utils.normalizeOptionalText(meta.halo)
        def t_min = Utils.normalizeOptionalText(meta.t_min)
        def t_max = Utils.normalizeOptionalText(meta.t_max)

        [
            id          : meta.id,
            zarr_output : meta.zarr_output,
            halo        : halo ? halo as float : null,
            t_min       : t_min ? t_min as int : null,
            t_max       : t_max ? t_max as int : null,
        ]
    }

    TILE_CONSTRUCTOR(input_ch.map { meta, dataset -> tuple(tileConstructorMeta(meta), dataset) })

    ch_ilastik_inputs = TILE_CONSTRUCTOR.out.tiles
        .map { meta, tiles_csv -> tuple(meta.id, tiles_csv) }
        .join(input_ch.map { meta, dataset -> tuple(meta.id, meta, dataset) })
        .flatMap { _id, tiles_csv, meta, raw_zarr ->
            def base_meta = ilastikMeta(meta)

            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = base_meta + [
                    resolution_level : row.resolution_level as Integer,
                    y_min : row.y_min as Integer,
                    y_max : row.y_max as Integer,
                    x_min : row.x_min as Integer,
                    x_max : row.x_max as Integer,
                    z_min : (row.z_min && row.z_min != '') ? row.z_min as Integer : null,
                    z_max : (row.z_max && row.z_max != '') ? row.z_max as Integer : null,
                ]

                return [ tile_meta, raw_zarr, meta.zarr_output ]
            }
        }
    ch_project = Channel.value(file(params.ilastik_project))

    ILASTIK_ZARR_PIPELINE(ch_ilastik_inputs, ch_project)
    ILASTIK_ZARR_PIPELINE.out.probabilities.view()
}
