#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP_ZARR                } from '../../modules/local/pymif_mip_zarr/main'
include { PYMIF_LABEL_MASK_ZARR         } from '../../modules/local/pymif_label_mask_zarr/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {

    def repoRootDir = workflow.projectDir.getParent().getParent()

    def default_zarr_dir = "${file(params.default_outdir)}/zarr"
    def default_mip_dir = "${file(params.default_outdir)}/mip"
    def default_mask_dir = "${file(params.default_outdir)}/masked"
    def default_labels_dir = "${file(params.default_outdir)}/labels"

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${workflow.projectDir}/assets/NO_MODEL_FILE")
    )

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    input_csv_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)

    pymif_inputs_ch = input_csv_ch
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            def requested_zarr_output = Utils.normalizeOptionalText(row.output)
            if (requested_zarr_output && !file(requested_zarr_output).isAbsolute()) {
                error(
                    "zarr_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_zarr_output}'. Leave the column blank to use the " +
                    "default location under default_outdir instead."
                )
            }
            def resolved_zarr_output = requested_zarr_output
                ?: "${default_zarr_dir}/${dataset_id}.ome.zarr"

            def requested_mip_output = Utils.normalizeOptionalText(row.mip_output)
            if (requested_mip_output && !file(requested_mip_output).isAbsolute()) {
                error(
                    "mip_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_mip_output}'. Leave the column blank to use the " +
                    "default location under default_outdir instead."
                )
            }
            def resolved_mip_output = requested_mip_output
                ?: "${default_mip_dir}/${dataset_id}_mip.ome.zarr"

            def requested_mask_output = Utils.normalizeOptionalText(row.mask_output)
            if (requested_mask_output && !file(requested_mask_output).isAbsolute()) {
                error(
                    "mask_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_mask_output}'. Leave the column blank to use the " +
                    "default location under default_outdir instead."
                )
            }
            def resolved_mask_output = requested_mask_output
                ?: "${default_mask_dir}/${dataset_id}_masked.ome.zarr"

            def requested_labels_output = Utils.normalizeOptionalText(row.labels_output)
            if (requested_labels_output && !file(requested_labels_output).isAbsolute()) {
                error(
                    "labels_output must be an absolute path (dataset_id=${dataset_id}): " +
                    "got '${requested_labels_output}'. Leave the column blank to use the " +
                    "default location under default_outdir instead."
                )
            }
            def resolved_labels_output = requested_labels_output
                ?: "${default_labels_dir}/${dataset_id}_labels.ome.zarr"

            def mask_channel_text = Utils.requireColumnValue(row, 'mask_channel')
            def mask_resolution_level = Utils.normalizeOptionalText(row.mask_resolution_level)
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)

            def diameter_value = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo = Utils.normalizeOptionalText(row.cellpose_halo)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)
            def do_3D_text = Utils.normalizeOptionalText(row.do_3D)
            def do_3D_value = do_3D_text ?: false

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                mip_output_name          : resolved_mip_output,

                mask_output_name         : resolved_mask_output,
                mask_channel             : mask_channel_text as int,
                mask_resolution_level    : mask_resolution_level
                    ? mask_resolution_level as int
                    : (resolution_level
                        ? resolution_level as int
                        : (params.resolution_level != null ? params.resolution_level as int : 0)),
                mask_min_area            : Utils.parseOptionalFloat(row.mask_min_area, 'mask_min_area'),
                mask_max_area            : Utils.parseOptionalFloat(row.mask_max_area, 'mask_max_area'),
                mask_intensity_channels  : Utils.normalizeOptionalText(row.mask_intensity_channels),
                mask_filters             : Utils.normalizeOptionalText(row.mask_filters),
                mask_otsu_multiplier     : Utils.normalizeOptionalText(row.mask_otsu_multiplier)
                    ? row.mask_otsu_multiplier as float
                    : 1.0,

                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                halo                  : cellpose_halo != null && cellpose_halo != ''
                    ? cellpose_halo as float
                    : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                use_physical_units    : use_physical_units ? use_physical_units.toBoolean() : false,
                do_3D                 : do_3D_value,
                labels_output_name    : resolved_labels_output
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    /*
     * Creates a multiscale 2D OME-Zarr by applying max() over Z at every
     * resolution level.
     *
     * Expected output:
     * tuple(meta, mip_zarr)
     */
    PYMIF_MIP_ZARR(PYMIF_CONVERSION.out.zarr)

    /*
     * Otsu-thresholds meta.mask_channel, labels + filters the resulting
     * objects (area and any other regionprops property via meta.mask_filters),
     * masks that channel to zero outside the kept objects, and writes a new
     * single-resolution OME-Zarr containing every original channel plus the
     * masked channel appended at the end.
     *
     * Expected output:
     * tuple(meta, masked_zarr)
     */
    PYMIF_LABEL_MASK_ZARR(PYMIF_MIP_ZARR.out.zarr)

    ch_cellpose_inputs = PYMIF_LABEL_MASK_ZARR.out.zarr
        .map { meta, zarr_path ->
            // PYMIF_LABEL_MASK_ZARR always writes a single-resolution OME-Zarr,
            // regardless of which mask_resolution_level/resolution_level the
            // CSV row requested upstream, so Cellpose must read level 0 of it.
            def cellpose_meta = meta + [resolution_level: 0]
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${workflow.projectDir}/assets/NO_FILE")
            )
            tuple(cellpose_meta, local_image, image_uri)
        }

    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)

}
