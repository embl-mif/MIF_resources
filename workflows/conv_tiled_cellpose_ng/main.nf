#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION                   } from '../../modules/local/pymif_conversion/main'
include { TILE_CONSTRUCTOR                   } from '../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF              } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {
    def zarr_dir = "${file(params.outdir)}/zarr"
 
    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->

            // Parameters for pymif conversion

            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = row.dataset_id?.toString()?.trim()
                ?: input_dataset.simpleName

            // Parameters for tile constructor

            def tile_size_x      = row.tile_size_x?.toString()?.trim()
            def tile_size_y      = row.tile_size_y?.toString()?.trim()
            def tile_size_z      = row.tile_size_z?.toString()?.trim()
            def tile_overlap     = row.tile_overlap?.toString()?.trim()
            def resolution_level = row.resolution_level?.toString()?.trim()
            def x_min            = row.x_min?.toString()?.trim()
            def x_max            = row.x_max?.toString()?.trim()
            def y_min            = row.y_min?.toString()?.trim()
            def y_max            = row.y_max?.toString()?.trim()
            def z_min            = row.z_min?.toString()?.trim()
            def z_max            = row.z_max?.toString()?.trim()
            def t_min            = row.t_min?.toString()?.trim()
            def t_max            = row.t_max?.toString()?.trim()
            def cellpose_halo    = params.cellpose_halo ?: 0.0

            // Parameters for cellpose segmentation

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

            // Construct a metadata map for the dataset, including all relevant parameters

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                tile_size_x      : tile_size_x ? tile_size_x as int : null,
                tile_size_y      : tile_size_y ? tile_size_y as int : null,
                tile_size_z      : tile_size_z ? tile_size_z as int : null,
                tile_overlap     : tile_overlap ? tile_overlap as float : null,
                resolution_level : resolution_level ? resolution_level as int : null,
                x_min            : x_min ? x_min as float : null,
                x_max            : x_max ? x_max as float : null,
                y_min            : y_min ? y_min as float : null,
                y_max            : y_max ? y_max as float : null,
                z_min            : z_min ? z_min as float : null,
                z_max            : z_max ? z_max as float : null,
                t_min            : t_min != null && t_min != '' ? t_min as int : null,
                t_max            : t_max != null && t_max != '' ? t_max as int : null,
                use_physical_units: false,
                cellpose_halo   : cellpose_halo,

                labels_output_name       : resolved_labels_output,
                diameter                 : diameter_value ? diameter_value as int : null,
                niter                    : niter_value ? niter_value as int : null,
                segmentation_channels    : cellpose_channels ?: null
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    TILE_CONSTRUCTOR(PYMIF_CONVERSION.out.zarr)

    ch_cellpose_inputs = TILE_CONSTRUCTOR.out.tiles
        .join(pymif_inputs_ch)
        .flatMap { meta, tiles_csv, raw_zarr ->
            tiles_csv.splitCsv(header: true).collect { row ->
                def tile_meta = meta.clone()
                tile_meta.y_min = row.y_min as Integer
                tile_meta.y_max = row.y_max as Integer
                tile_meta.x_min = row.x_min as Integer
                tile_meta.x_max = row.x_max as Integer
                tile_meta.z_min = (row.z_min && row.z_min != '') ? row.z_min as Integer : null
                tile_meta.z_max = (row.z_max && row.z_max != '') ? row.z_max as Integer : null
                tile_meta.halo = row.cellpose_halo
                return [ tile_meta, raw_zarr]
            }
        }
    CELLPOSE_SAM_ZARR(ch_cellpose_inputs)
}

