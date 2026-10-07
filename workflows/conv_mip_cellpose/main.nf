#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION      } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP             } from '../../modules/local/pymif_mip/main'
include { CELLPOSE_SAM_ZARR_MIF } from '../../modules/local/cellpose_sam_zarr_mif/main'


workflow {

    def repoRootDir = workflow.projectDir.getParent().getParent()

    // Only evaluated for a row that leaves an output column blank.
    def defaultOutput = { String subdir, String fileName ->
        if (!params.outdir) {
            error "--outdir is required when a row leaves its output column blank (needed for ${fileName})"
        }
        "${file(params.outdir)}/${subdir}/${fileName}".toString()
    }

    // Absolute path from the row if given, otherwise <outdir>/<subdir>/<fileName>.
    def resolveRowOutput = { row, String column, String subdir, String fileName, String datasetId ->
        Utils.normalizeOptionalText(row[column])
            ? Utils.resolveOutput(row, column, null, datasetId)
            : defaultOutput(subdir, fileName)
    }

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${repoRootDir}/assets/NO_MODEL_FILE")
    )

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    // Stage 1: conversion. All three output paths are resolved here so a bad or
    // missing path fails before any task runs. Each stage then overrides
    // zarr_output on the meta it consumes, so no stage writes over the previous one.
    pymif_inputs_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def zarr_output   = resolveRowOutput(row, 'output',        'zarr',   "${dataset_id}.ome.zarr",        dataset_id)
            def mip_output    = resolveRowOutput(row, 'mip_output',    'mip',    "${dataset_id}_mip.ome.zarr",    dataset_id)
            def labels_output = resolveRowOutput(row, 'labels_output', 'labels', "${dataset_id}_labels.ome.zarr", dataset_id)

            def meta = row + [
                id           : dataset_id,
                zarr_output  : zarr_output,
                mip_output   : mip_output,
                labels_output: labels_output,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    // Stage 2: MIP writes its own store.
    PYMIF_MIP(
        PYMIF_CONVERSION.out.zarr.map { meta, zarr_path ->
            tuple(meta + [zarr_output: meta.mip_output], zarr_path)
        }
    )

    // Stage 3: Cellpose. Own meta with only the keys the module reads, and
    // zarr_output set to the labels path.
    ch_cellpose_inputs = PYMIF_MIP.out.zarr
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${repoRootDir}/assets/NO_FILE")
            )

            def resolution_level = Utils.parseOptionalInt(meta.resolution_level, 'resolution_level')
            def halo             = Utils.parseOptionalFloat(meta.cellpose_halo, 'cellpose_halo')
            def physical         = Utils.normalizeOptionalText(meta.use_physical_units)
            def do3d             = Utils.normalizeOptionalText(meta.do_3D)

            def cellposeMeta = [
                id                    : meta.id,
                zarr_output           : meta.labels_output,
                resolution_level      : resolution_level != null ? resolution_level : (params.resolution_level != null ? params.resolution_level as int : null),
                use_physical_units    : physical ? physical.toBoolean() : false,
                do_3D                 : do3d ? do3d.equalsIgnoreCase('true') : false,
                diameter              : Utils.parseOptionalInt(meta.cellpose_diameter, 'cellpose_diameter'),
                niter                 : Utils.parseOptionalInt(meta.cellpose_niter, 'cellpose_niter'),
                segmentation_channels : Utils.normalizeOptionalText(meta.cellpose_channels),
                flow_threshold        : Utils.parseOptionalFloat(meta.cellpose_flow_threshold, 'cellpose_flow_threshold'),
                cellprob_threshold    : Utils.parseOptionalFloat(meta.cellpose_cellprob_threshold, 'cellpose_cellprob_threshold'),
                halo                  : halo != null ? halo : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                x_min                 : Utils.parseOptionalFloat(meta.x_min, 'x_min'),
                x_max                 : Utils.parseOptionalFloat(meta.x_max, 'x_max'),
                y_min                 : Utils.parseOptionalFloat(meta.y_min, 'y_min'),
                y_max                 : Utils.parseOptionalFloat(meta.y_max, 'y_max'),
                z_min                 : Utils.parseOptionalFloat(meta.z_min, 'z_min'),
                z_max                 : Utils.parseOptionalFloat(meta.z_max, 'z_max'),
            ]

            tuple(cellposeMeta, local_image, image_uri)
        }

    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)
}
