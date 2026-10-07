#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION      } from '../../modules/local/pymif_conversion/main'
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

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${repoRootDir}/assets/NO_MODEL_FILE")
    )

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    // Stage 1: conversion. Resolves both output paths up front so a bad or missing
    // path fails before any task runs. The row keeps its raw columns; each later
    // stage parses only the ones it consumes.
    pymif_inputs_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def zarr_output = Utils.normalizeOptionalText(row.output)
                ? Utils.resolveOutput(row, 'output', null, dataset_id)
                : defaultOutput('zarr', "${dataset_id}.ome.zarr")

            def labels_output = Utils.normalizeOptionalText(row.labels_output)
                ? Utils.resolveOutput(row, 'labels_output', null, dataset_id)
                : defaultOutput('labels', "${dataset_id}_labels.ome.zarr")

            def meta = row + [
                id           : dataset_id,
                zarr_output  : zarr_output,
                labels_output: labels_output,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    // Stage 2: Cellpose. Own meta with only the keys the module reads, and its own
    // zarr_output so it never writes over the converted image.
    ch_cellpose_inputs = PYMIF_CONVERSION.out.zarr
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${repoRootDir}/assets/NO_FILE")
            )

            def intOrParam = { column, fallback ->
                def v = Utils.parseOptionalInt(meta[column], column)
                v != null ? v : (fallback != null ? fallback as int : null)
            }
            def halo = Utils.parseOptionalFloat(meta.cellpose_halo, 'cellpose_halo')
            def physical = Utils.normalizeOptionalText(meta.use_physical_units)
            def do3d = Utils.normalizeOptionalText(meta.do_3D)

            def cellposeMeta = [
                id                    : meta.id,
                zarr_output           : meta.labels_output,
                resolution_level      : intOrParam('resolution_level', params.resolution_level),
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
