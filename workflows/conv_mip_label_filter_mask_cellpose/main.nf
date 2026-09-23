#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP                     } from '../../modules/local/pymif_mip/main'
include { PYMIF_OTSU_LABEL              } from '../../modules/local/pymif_otsu_label/main'
include { PYMIF_FEATURE_EXTRACTION      } from '../../modules/local/pymif_feature_extraction/main'
include { PYMIF_FILTER                  } from '../../modules/local/pymif_filter/main'
include { PYMIF_CLEAN_CHANNEL           } from '../../modules/local/pymif_clean_channel/main'
include { CELLPOSE_SAM_ZARR_MIF         } from '../../modules/local/cellpose_sam_zarr_mif/main'

workflow {

    def default_zarr_dir = "${file(params.default_outdir)}/zarr"
    def default_mip_dir = "${file(params.default_outdir)}/mip"
    def default_otsu_dir = "${file(params.default_outdir)}/otsu_labels"
    def default_features_dir = "${file(params.default_outdir)}/measurements"
    def default_filtered_dir = "${file(params.default_outdir)}/filtered"
    def default_clean_dir = "${file(params.default_outdir)}/clean"
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

            def resolveOutput = { String columnName, String defaultPath ->
                def requested = Utils.normalizeOptionalText(row[columnName])
                if (requested && !file(requested).isAbsolute()) {
                    error(
                        "${columnName} must be an absolute path (dataset_id=${dataset_id}): " +
                        "got '${requested}'. Leave the column blank to use the default " +
                        "location under default_outdir instead."
                    )
                }
                requested ?: defaultPath
            }

            def resolved_zarr_output = resolveOutput('output', "${default_zarr_dir}/${dataset_id}.ome.zarr")
            def resolved_mip_output = resolveOutput('mip_output', "${default_mip_dir}/${dataset_id}_mip.ome.zarr")
            def resolved_otsu_output = resolveOutput('otsu_output', "${default_otsu_dir}/${dataset_id}_otsu_labels.ome.zarr")
            def resolved_features_output = resolveOutput('features_output', "${default_features_dir}/${dataset_id}_features.csv")
            def resolved_filter_output = resolveOutput('filter_output', "${default_filtered_dir}/${dataset_id}_filtered_labels.ome.zarr")
            def resolved_filter_csv_output = resolveOutput('filter_csv_output', "${default_filtered_dir}/${dataset_id}_filtered_labels.csv")
            def resolved_clean_output = resolveOutput('clean_output', "${default_clean_dir}/${dataset_id}_clean.ome.zarr")
            def resolved_labels_output = resolveOutput('labels_output', "${default_labels_dir}/${dataset_id}_labels.ome.zarr")

            def otsu_channel_text = Utils.requireColumnValue(row, 'otsu_channel')
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)
            def otsu_resolution_level_text = Utils.normalizeOptionalText(row.otsu_resolution_level)
            def resolved_otsu_resolution_level = otsu_resolution_level_text
                ? otsu_resolution_level_text as int
                : (resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : 0))

            def features_resolution_level_text = Utils.normalizeOptionalText(row.features_resolution_level)
            def resolved_features_resolution_level = features_resolution_level_text
                ? features_resolution_level_text as int
                : resolved_otsu_resolution_level

            def clean_channel_text = Utils.normalizeOptionalText(row.clean_channel)

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

                // PYMIF_OTSU_LABEL: threshold + label meta.otsu_channel of the
                // MIP, rebuilding the full pyramid from otsu_resolution_level.
                otsu_output_name         : resolved_otsu_output,
                otsu_channel             : otsu_channel_text as int,
                otsu_resolution_level    : resolved_otsu_resolution_level,
                otsu_multiplier          : Utils.normalizeOptionalText(row.otsu_multiplier)
                    ? row.otsu_multiplier as float
                    : 1.0,

                // PYMIF_FEATURE_EXTRACTION: measure the MIP + Otsu labels.
                // Both zarrs share the same pyramid, so the same level is
                // used for both unless features_resolution_level overrides it.
                features_output_name                : resolved_features_output,
                features_intensity_resolution_level : resolved_features_resolution_level,
                features_labels_resolution_level    : resolved_features_resolution_level,
                features_channels                   : Utils.normalizeOptionalText(row.features_channels),
                features_extra_properties           : Utils.normalizeOptionalText(row.features_extra_properties),
                features_extra_intensity_properties : Utils.normalizeOptionalText(row.features_extra_intensity_properties),

                // PYMIF_FILTER: keep only labels whose measured row satisfies
                // these limits, and relabel the survivors to 1..N.
                filter_output_name      : resolved_filter_output,
                filter_csv_output_name  : resolved_filter_csv_output,
                filter_min_area         : Utils.parseOptionalFloat(row.filter_min_area, 'filter_min_area'),
                filter_max_area         : Utils.parseOptionalFloat(row.filter_max_area, 'filter_max_area'),
                filter_properties       : Utils.normalizeOptionalText(row.filter_properties),

                // PYMIF_CLEAN_CHANNEL: zero meta.clean_channel of the MIP
                // outside the filtered labels, and append it as a new channel.
                clean_output_name        : resolved_clean_output,
                clean_channel            : clean_channel_text
                    ? clean_channel_text as int
                    : otsu_channel_text as int,

                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                halo                  : cellpose_halo != null && cellpose_halo != ''
                    ? cellpose_halo as float
                    : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                resolution_level      : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : 0),
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
    PYMIF_MIP(PYMIF_CONVERSION.out.zarr)

    /*
     * Otsu-thresholds + labels meta.otsu_channel of the MIP at
     * meta.otsu_resolution_level, then rebuilds the full pyramid (same
     * num_levels/chunksize/downsample factors as the MIP) by nearest-neighbor
     * resizing that canonical labeling to every other level.
     *
     * Expected output:
     * tuple(meta, otsu_labels_zarr)
     */
    PYMIF_OTSU_LABEL(PYMIF_MIP.out.zarr)

    /*
     * Measures every Otsu-labeled object against the MIP (every channel by
     * default) with skimage.measure.regionprops.
     *
     * Expected output:
     * tuple(meta, features_csv)
     */
    ch_features_inputs = PYMIF_MIP.out.zarr
        .map { meta, zarr -> tuple(meta.id, meta, zarr) }
        .join(PYMIF_OTSU_LABEL.out.zarr.map { meta, labels -> tuple(meta.id, labels) })
        .map { id, meta, intensity, labels -> tuple(meta, intensity, labels) }

    PYMIF_FEATURE_EXTRACTION(ch_features_inputs)

    /*
     * Keeps only the Otsu-labeled objects whose measured row satisfies
     * meta.filter_min_area/filter_max_area/filter_properties, relabels the
     * survivors to a contiguous 1..N id space, and writes both the filtered
     * labels and the filtered measurements CSV.
     *
     * Expected output:
     * tuple(meta, filtered_labels_zarr, filtered_csv)
     */
    ch_filter_inputs = PYMIF_OTSU_LABEL.out.zarr
        .map { meta, labels -> tuple(meta.id, meta, labels) }
        .join(PYMIF_FEATURE_EXTRACTION.out.csv.map { meta, csv -> tuple(meta.id, csv) })
        .map { id, meta, labels, csv -> tuple(meta, labels, csv) }

    PYMIF_FILTER(ch_filter_inputs)

    /*
     * Zeros meta.clean_channel of the MIP outside the filtered labels, and
     * appends it as a new channel of a new OME-Zarr (every original MIP
     * channel plus this one), preserving the MIP's full pyramid.
     *
     * Expected output:
     * tuple(meta, clean_zarr)
     */
    ch_clean_inputs = PYMIF_MIP.out.zarr
        .map { meta, intensity -> tuple(meta.id, meta, intensity) }
        .join(PYMIF_FILTER.out.zarr_csv.map { meta, labels, csv -> tuple(meta.id, labels) })
        .map { id, meta, intensity, labels -> tuple(meta, intensity, labels) }

    PYMIF_CLEAN_CHANNEL(ch_clean_inputs)

    ch_cellpose_inputs = PYMIF_CLEAN_CHANNEL.out.zarr
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${workflow.projectDir}/assets/NO_FILE")
            )
            tuple(meta, local_image, image_uri)
        }

    CELLPOSE_SAM_ZARR_MIF(ch_cellpose_inputs, ch_cellpose_model)

}
