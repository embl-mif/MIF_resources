#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { TILE_CONSTRUCTOR              } from '../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {
    def zarr_dir = "${file(params.outdir)}/zarr"
    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${workflow.projectDir}/assets/NO_MODEL_FILE")
    )

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            // Read the input dataset and prepare metadata for processing
            def input_dataset = file(row.input, checkIfExists: true)
            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def tile_size_x = Utils.normalizeOptionalText(row.tile_size_x)
            def tile_size_y = Utils.normalizeOptionalText(row.tile_size_y)
            def tile_size_z = Utils.normalizeOptionalText(row.tile_size_z)
            def tile_overlap = Utils.normalizeOptionalText(row.tile_overlap)
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)
            def x_min = Utils.normalizeOptionalText(row.x_min)
            def x_max = Utils.normalizeOptionalText(row.x_max)
            def y_min = Utils.normalizeOptionalText(row.y_min)
            def y_max = Utils.normalizeOptionalText(row.y_max)
            def z_min = Utils.normalizeOptionalText(row.z_min)
            def z_max = Utils.normalizeOptionalText(row.z_max)
            def t_min = Utils.normalizeOptionalText(row.t_min)
            def t_max = Utils.normalizeOptionalText(row.t_max)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)

            def diameter_value = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo = Utils.normalizeOptionalText(row.cellpose_halo)

            def do_3D_text = Utils.normalizeOptionalText(row.do_3D)
            def do_3D_value = do_3D_text
                ? do_3D_text.equalsIgnoreCase('true')
                : false

            def resolved_zarr_output =
                "${zarr_dir}/${file(row.output).name}"
            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                tile_size_x      : tile_size_x ? tile_size_x as int : null,
                tile_size_y      : tile_size_y ? tile_size_y as int : null,
                tile_size_z      : tile_size_z ? tile_size_z as int : null,
                tile_overlap     : tile_overlap ? tile_overlap as float : null,
                resolution_level      : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : null),
                x_min            : x_min ? x_min as float : null,
                x_max            : x_max ? x_max as float : null,
                y_min            : y_min ? y_min as float : null,
                y_max            : y_max ? y_max as float : null,
                z_min            : z_min ? z_min as float : null,
                z_max            : z_max ? z_max as float : null,
                t_min            : t_min != null && t_min != '' ? t_min as int : null,
                t_max            : t_max != null && t_max != '' ? t_max as int : null,
                use_physical_units: use_physical_units ? use_physical_units.toBoolean() : false,
                diameter                 : diameter_value ? diameter_value as int : null,
                niter                    : niter_value ? niter_value as int : null,
                segmentation_channels    : cellpose_channels ?: null,
                halo                     : cellpose_halo != null && cellpose_halo != '' ? cellpose_halo as float : 20.0,
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
                ]
                def (local_image, image_uri) = Utils.splitImageInput(
                    meta.zarr_output_name,
                    file("${workflow.projectDir}/assets/NO_FILE")
                )
                tuple(tile_meta, local_image, image_uri)
            }
        }
    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)
}

