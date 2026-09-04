#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { TILECONSTRUCTOR_INPUT_PREP    } from '../../modules/local/tileconstructor_input_prep/main'
include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {
    def zarr_dir = "${file(params.outdir)}/zarr"
    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->

            def index = row_index.getAndIncrement()

            // Read the input dataset and prepare metadata for processing
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = row.dataset_id?.toString()?.trim()
                ?: "${input_dataset.simpleName}_${index + 1}"

            // Read optional parameters for tiling
            def tile_size_x = row.tile_size_x?.toString()?.trim()
            def tile_size_y = row.tile_size_y?.toString()?.trim()
            def tile_size_z = row.tile_size_z?.toString()?.trim()
            def tile_overlap = row.tile_overlap?.toString()?.trim()
            def resolution_level = row.resolution_level?.toString()?.trim()
            def x_min = row.x_min?.toString()?.trim()
            def x_max = row.x_max?.toString()?.trim()
            def y_min = row.y_min?.toString()?.trim()
            def y_max = row.y_max?.toString()?.trim()
            def z_min = row.z_min?.toString()?.trim()
            def z_max = row.z_max?.toString()?.trim()
            def t_min = row.t_min?.toString()?.trim()
            def t_max = row.t_max?.toString()?.trim()
            def use_physical_units = row.use_physical_units?.toString()?.trim()

            // Read optional parameters for Cellpose processing
            def diameter_value = row.cellpose_diameter?.toString()?.trim()
            def niter_value = row.cellpose_niter?.toString()?.trim()
            def cellpose_channels = row.cellpose_channels?.toString()?.trim()
            def cellpose_halo = row.cellpose_halo?.toString()?.trim()

            def do_3D_text = row.do_3D?.toString()?.trim()
            def do_3D_value = do_3D_text
                ? do_3D_text.equalsIgnoreCase('true')
                : false

            def resolved_zarr_output =
                "${zarr_dir}/${file(row.output).name}"
            
            def labels_output_value = row.labels_output?.toString()?.trim()
            def labels_output_name = labels_output_value
                ? file(labels_output_value).name
                : "${dataset_id}_labels.ome.zarr"
            def resolved_labels_output = "${zarr_dir}/${labels_output_name}"

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                tile_size_x      : tile_size_x ? tile_size_x as int : 512,
                tile_size_y      : tile_size_y ? tile_size_y as int : 512,
                tile_size_z      : tile_size_z ? tile_size_z as int : 512,
                tile_overlap     : tile_overlap ? tile_overlap as float : 32,
                resolution_level : resolution_level ? resolution_level as int : null,
                x_min            : x_min ? x_min as float : null,
                x_max            : x_max ? x_max as float : null,
                y_min            : y_min ? y_min as float : null,
                y_max            : y_max ? y_max as float : null,
                z_min            : z_min ? z_min as float : null,
                z_max            : z_max ? z_max as float : null,
                t_min            : t_min != null && t_min != '' ? t_min as int : null,
                t_max            : t_max != null && t_max != '' ? t_max as int : null,
                use_physical_units: use_physical_units ? use_physical_units.toBoolean() : false,

                labels_output_name       : resolved_labels_output,
                diameter                 : diameter_value ? diameter_value as int : null,
                niter                    : niter_value ? niter_value as int : null,
                segmentation_channels    : cellpose_channels ?: null,
                cellpose_halo            : cellpose_halo != null && cellpose_halo != '' ? cellpose_halo as float : 20.0,
                do_3D                    : do_3D_value
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    TILE_CONSTRUCTOR(PYMIF_CONVERSION.out.zarr)

    ch_cellpose_inputs = TILE_CONSTRUCTOR.out.tiles
        .flatMap { meta, tiles_csv ->
            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta + [
                    id                 : row.dataset_id,
                    dataset_id         : row.dataset_id,
                    source_dataset_id  : row.source_dataset_id,
                    tile_id            : row.tile_id,
                    tile_index         : row.tile_index as Integer,

                    x_min              : row.x_min as Integer,
                    x_max              : row.x_max as Integer,
                    y_min              : row.y_min as Integer,
                    y_max              : row.y_max as Integer,
                    z_min              : row.z_min ? row.z_min as Integer : null,
                    z_max              : row.z_max ? row.z_max as Integer : null,

                    labels_output_name : "${zarr_dir}/${row.dataset_id}_labels.ome.zarr"
                ]

                tuple(tile_meta, meta.zarr_output_name)
            }
        }
    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs)
}

