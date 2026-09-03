#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { TILECONSTRUCTOR_INPUT_PREP    } from '../../modules/local/tileconstructor_input_prep/main'
include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {
    def zarr_dir = "${file(params.outdir)}/zarr"

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = row.dataset_id?.toString()?.trim()
                ?: input_dataset.simpleName

            def diameter_value = row.cellpose_diameter?.toString()?.trim()
            def niter_value = row.cellpose_niter?.toString()?.trim()
            def cellpose_channels = row.cellpose_channels?.toString()?.trim()

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
                labels_output_name       : resolved_labels_output,
                diameter                 : diameter_value ? diameter_value as int : null,
                niter                    : niter_value ? niter_value as int : null,
                segmentation_channels    : cellpose_channels ?: null
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    ch_zarr_conversion_done = PYMIF_CONVERSION.out.zarr.collect()

    config_csv_ch = TILECONSTRUCTOR_INPUT_PREP(ch_zarr_conversion_done.map { zarr_dir }, params.resolution_level)

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
    def cellpose_halo = params.cellpose_halo ?: 0.0

    ch_cellpose_input = input_ch.map { meta, image -> tuple(meta + [halo: cellpose_halo], image) }

    TILE_CONSTRUCTOR(ch_cellpose_input)

    ch_cellpose_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(ch_cellpose_input)
        .flatMap { meta, tiles_csv, raw_zarr ->
            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta.clone()
                tile_meta.y_min = row.y_min as Integer
                tile_meta.y_max = row.y_max as Integer
                tile_meta.x_min = row.x_min as Integer
                tile_meta.x_max = row.x_max as Integer
                tile_meta.z_min = (row.z_min && row.z_min != '') ? row.z_min as Integer : null
                tile_meta.z_max = (row.z_max && row.z_max != '') ? row.z_max as Integer : null
                tile_meta.halo = cellpose_halo
                return [ tile_meta, raw_zarr]
            }
        }
    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs)
}

