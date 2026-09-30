#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION              } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP                     } from '../../modules/local/pymif_mip/main'
include { PYMIF_FEATURE_EXTRACTION      } from '../../modules/local/pymif_feature_extraction/main'
include { PYMIF_FILTER                  } from '../../modules/local/pymif_filter/main'
include { PYMIF_CLEAN_CHANNEL           } from '../../modules/local/pymif_clean_channel/main'
include { CELLPOSE_SAM_ZARR_MIF as CELLPOSE_SEED  } from '../../modules/local/cellpose_sam_zarr_mif/main'
include { CELLPOSE_SAM_ZARR_MIF as CELLPOSE_FINAL } from '../../modules/local/cellpose_sam_zarr_mif/main'
include { PYMIF_ENDOMORPH_REGIONS       } from '../../modules/local/pymif_endomorph_regions/main'
include { PYMIF_ENDOMORPH_FEATURES      } from '../../modules/local/pymif_endomorph_features/main'
include { PYMIF_ENDOMORPH_MERGE         } from '../../modules/local/pymif_endomorph_merge/main'

workflow {

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

            def resolved_zarr_output = Utils.resolveOutput(row, 'output', "${dataset_id}.ome.zarr", dataset_id)
            def resolved_mip_output = Utils.resolveOutput(row, 'mip_output', "${dataset_id}_mip.ome.zarr", dataset_id)
            def resolved_seed_output = Utils.resolveOutput(row, 'seed_output', "${dataset_id}_seed_labels.ome.zarr", dataset_id)
            def resolved_features_output = Utils.resolveOutput(row, 'features_output', "${dataset_id}_features.csv", dataset_id)
            def resolved_filter_output = Utils.resolveOutput(row, 'filter_output', "${dataset_id}_filtered_labels.ome.zarr", dataset_id)
            def resolved_filter_csv_output = Utils.resolveOutput(row, 'filter_csv_output', "${dataset_id}_filtered_labels.csv", dataset_id)
            def resolved_clean_output = Utils.resolveOutput(row, 'clean_output', "${dataset_id}_clean.ome.zarr", dataset_id)
            def resolved_labels_output = Utils.resolveOutput(row, 'labels_output', "${dataset_id}_labels.ome.zarr", dataset_id)
            def resolved_endomorph_regions_output = Utils.resolveOutput(row, 'endomorph_regions_output', "${dataset_id}_endomorph_regions.ome.zarr", dataset_id)
            def resolved_endomorph_cells_output = Utils.resolveOutput(row, 'endomorph_cells_output', "${dataset_id}_endomorph_cells.csv", dataset_id)
            def resolved_endomorph_summary_output = Utils.resolveOutput(row, 'endomorph_summary_output', "${dataset_id}_endomorph_summary.csv", dataset_id)

            // Required: which MIP channel(s) the seed Cellpose call segments
            // on, in place of the Otsu-threshold channel used by
            // conv_mip_label_filter_mask_cellpose.
            def seed_channels_text = Utils.requireColumnValue(row, 'seed_channels')
            def resolution_level = Utils.normalizeOptionalText(row.resolution_level)

            def features_resolution_level_text = Utils.normalizeOptionalText(row.features_resolution_level)
            // The seed Cellpose call always segments at resolution level 0
            // (see the PYMIF_FEATURE_EXTRACTION note below), so that is also
            // the only resolution level guaranteed valid for both the MIP
            // and the seed labels here.
            def resolved_features_resolution_level = features_resolution_level_text
                ? features_resolution_level_text as int
                : 0

            def clean_channel_text = Utils.normalizeOptionalText(row.clean_channel)
            // No single default channel exists when seed_channels lists
            // several channels, so fall back to the first one listed.
            def default_clean_channel = (seed_channels_text.split(/[;,\s]+/).findAll { it })[0] as int

            def seed_diameter_value = Utils.normalizeOptionalText(row.seed_diameter)
            def seed_niter_value = Utils.normalizeOptionalText(row.seed_niter)
            def seed_halo = Utils.normalizeOptionalText(row.seed_halo)
            def seed_flow_threshold = Utils.normalizeOptionalText(row.seed_flow_threshold)
            def seed_cellprob_threshold = Utils.normalizeOptionalText(row.seed_cellprob_threshold)
            def seed_use_physical_units = Utils.normalizeOptionalText(row.seed_use_physical_units)
            def seed_do_3D_text = Utils.normalizeOptionalText(row.seed_do_3D)

            def diameter_value = Utils.normalizeOptionalText(row.cellpose_diameter)
            def niter_value = Utils.normalizeOptionalText(row.cellpose_niter)
            def cellpose_channels = Utils.normalizeOptionalText(row.cellpose_channels)
            def cellpose_halo = Utils.normalizeOptionalText(row.cellpose_halo)
            def flow_threshold = Utils.normalizeOptionalText(row.cellpose_flow_threshold)
            def cellprob_threshold = Utils.normalizeOptionalText(row.cellpose_cellprob_threshold)
            def use_physical_units = Utils.normalizeOptionalText(row.use_physical_units)
            def do_3D_text = Utils.normalizeOptionalText(row.do_3D)

            def endomorph_border_margin_value = Utils.normalizeOptionalText(row.endomorph_border_margin)
            def endomorph_timepoint_value = Utils.normalizeOptionalText(row.endomorph_timepoint)
            def endomorph_empty_smoothing_value = Utils.normalizeOptionalText(row.endomorph_empty_smoothing)
            def endomorph_empty_min_area_value = Utils.normalizeOptionalText(row.endomorph_empty_min_area)
            def endomorph_empty_seed_percentile_value = Utils.normalizeOptionalText(row.endomorph_empty_seed_percentile)
            def endomorph_empty_min_contrast_value = Utils.normalizeOptionalText(row.endomorph_empty_min_contrast)

            def meta = row + [
                id                       : dataset_id,
                zarr_output_name         : resolved_zarr_output,

                mip_output_name          : resolved_mip_output,

                // CELLPOSE_SEED: segments meta.seed_segmentation_channels of
                // the MIP at resolution level 0, replacing the Otsu-label
                // step. Kept under seed_*-prefixed keys (rather than the
                // plain segmentation_channels/halo/... keys CELLPOSE_SAM_ZARR_MIF
                // reads) so they don't collide with the final Cellpose call's
                // own settings below; they're remapped just before that call.
                seed_output_name           : resolved_seed_output,
                seed_segmentation_channels : seed_channels_text,
                seed_halo                  : seed_halo != null && seed_halo != '' ? seed_halo as float : 20.0,
                seed_diameter               : seed_diameter_value ? seed_diameter_value as int : null,
                seed_niter                  : seed_niter_value ? seed_niter_value as int : null,
                seed_flow_threshold         : seed_flow_threshold ? seed_flow_threshold as float : null,
                seed_cellprob_threshold     : seed_cellprob_threshold ? seed_cellprob_threshold as float : null,
                seed_use_physical_units     : seed_use_physical_units ? seed_use_physical_units.toBoolean() : false,
                seed_do_3D                  : seed_do_3D_text ?: false,

                // PYMIF_FEATURE_EXTRACTION: measure the MIP + seed Cellpose
                // labels. Both zarrs share the same pyramid only when the
                // seed call ran at resolution level 0 (CELLPOSE_SAM_ZARR_MIF
                // otherwise builds a shorter pyramid starting at that level),
                // so this is not exposed as a configurable column here.
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
                    : default_clean_channel,

                // CELLPOSE_FINAL: segments the clean-channel OME-Zarr.
                diameter              : diameter_value ? diameter_value as int : null,
                niter                 : niter_value ? niter_value as int : null,
                segmentation_channels : cellpose_channels ?: null,
                halo                  : cellpose_halo != null && cellpose_halo != ''
                    ? cellpose_halo as float
                    : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
                flow_threshold        : flow_threshold ? flow_threshold as float : null,
                cellprob_threshold    : cellprob_threshold ? cellprob_threshold as float : null,
                resolution_level      : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : 0),
                use_physical_units    : use_physical_units ? use_physical_units.toBoolean() : false,
                do_3D                 : do_3D_text ?: false,
                labels_output_name    : resolved_labels_output,

                // PYMIF_ENDOMORPH_REGIONS: for every final-Cellpose-labelled
                // cell that doesn't touch the image border, computes nested
                // R1/R2/R3 edge-distance bands (inward from the cell edge)
                // as a label OME-Zarr, and estimates the cell-free ("empty")
                // area from the MIP fluorescence signal (labels 4=empty,
                // 5=not imaged). Morphology and intensities are computed
                // downstream by PYMIF_ENDOMORPH_FEATURES.
                endomorph_regions_output_name     : resolved_endomorph_regions_output,
                endomorph_ring_depths             : Utils.normalizeOptionalText(row.endomorph_ring_depths) ?: '6 9 24',
                endomorph_border_margin           : endomorph_border_margin_value ? endomorph_border_margin_value as int : 10,
                endomorph_timepoint               : endomorph_timepoint_value ? endomorph_timepoint_value as int : 0,

                // PYMIF_ENDOMORPH_FEATURES: measures per-cell morphology and
                // V-cadherin/myosin/actin band intensities against the bands
                // above, and estimates the cell-free area of the image from
                // the fluorescence signal. CELLPOSE_FINAL's output labels
                // pyramid starts at meta.resolution_level (see the seed
                // features note above), so its own level 0 is the only one
                // guaranteed to share the MIP's shape at that level.
                endomorph_intensity_resolution_level : resolution_level
                    ? resolution_level as int
                    : (params.resolution_level != null ? params.resolution_level as int : 0),
                endomorph_labels_resolution_level    : 0,
                endomorph_vcad_channel        : Utils.normalizeOptionalText(row.endomorph_vcad_channel),
                endomorph_myosin_channel      : Utils.normalizeOptionalText(row.endomorph_myosin_channel),
                endomorph_actin_channel       : Utils.normalizeOptionalText(row.endomorph_actin_channel),
                endomorph_nuclei_channel      : Utils.normalizeOptionalText(row.endomorph_nuclei_channel),
                endomorph_empty_channels          : Utils.normalizeOptionalText(row.endomorph_empty_channels) ?: 'vcad myosin actin',
                endomorph_empty_smoothing         : endomorph_empty_smoothing_value ? endomorph_empty_smoothing_value as float : 2.0,
                endomorph_empty_min_area          : endomorph_empty_min_area_value ? endomorph_empty_min_area_value as float : 25.0,
                endomorph_empty_seed_percentile   : endomorph_empty_seed_percentile_value ? endomorph_empty_seed_percentile_value as float : 1.0,
                endomorph_empty_min_contrast      : endomorph_empty_min_contrast_value ? endomorph_empty_min_contrast_value as float : 0.25,
                endomorph_cells_output_name       : resolved_endomorph_cells_output,
                endomorph_summary_output_name     : resolved_endomorph_summary_output
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
     * Segments meta.seed_segmentation_channels of the MIP with Cellpose-SAM
     * at resolution level 0, in place of the Otsu-threshold + label step.
     * Remap the seed_*-prefixed meta fields into the plain field names
     * CELLPOSE_SAM_ZARR_MIF reads, and force resolution_level 0 so its
     * output pyramid has exactly as many levels as the MIP (see the
     * features_resolution_level note above).
     *
     * Expected output:
     * tuple(meta, seed_labels_zarr)
     */
    ch_seed_inputs = PYMIF_MIP.out.zarr
        .map { meta, zarr_path ->
            def seed_meta = meta + [
                segmentation_channels: meta.seed_segmentation_channels,
                halo                 : meta.seed_halo,
                diameter             : meta.seed_diameter,
                niter                : meta.seed_niter,
                flow_threshold       : meta.seed_flow_threshold,
                cellprob_threshold   : meta.seed_cellprob_threshold,
                use_physical_units   : meta.seed_use_physical_units,
                do_3D                : meta.seed_do_3D,
                resolution_level     : 0,
                labels_output_name   : meta.seed_output_name,
            ]
            def (local_image, image_uri) = Utils.splitImageInput(
                zarr_path,
                file("${workflow.projectDir}/assets/NO_FILE")
            )
            tuple(seed_meta, local_image, image_uri)
        }

    CELLPOSE_SEED(ch_seed_inputs, ch_cellpose_model)

    /*
     * Measures every seed-labeled object against the MIP (every channel by
     * default) with skimage.measure.regionprops.
     *
     * Expected output:
     * tuple(meta, features_csv)
     */
    ch_features_inputs = PYMIF_MIP.out.zarr
        .map { meta, zarr -> tuple(meta.id, meta, zarr) }
        .join(CELLPOSE_SEED.out.labels.map { meta, labels -> tuple(meta.id, labels) })
        .map { id, meta, intensity, labels -> tuple(meta, intensity, labels) }

    PYMIF_FEATURE_EXTRACTION(ch_features_inputs)

    /*
     * Keeps only the seed-labeled objects whose measured row satisfies
     * meta.filter_min_area/filter_max_area/filter_properties, relabels the
     * survivors to a contiguous 1..N id space, and writes both the filtered
     * labels and the filtered measurements CSV.
     *
     * Expected output:
     * tuple(meta, filtered_labels_zarr, filtered_csv)
     */
    ch_filter_inputs = CELLPOSE_SEED.out.labels
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

    CELLPOSE_FINAL(ch_cellpose_inputs, ch_cellpose_model)

    /*
     * For every final-Cellpose-labelled cell that doesn't touch the image
     * border, computes nested R1/R2/R3 edge-distance bands (inward from the
     * cell edge) plus the empty (4) and not-imaged (5) areas of the MIP, as a
     * label OME-Zarr (0/1/2/3/4/5). A cell with no painted band pixel was
     * excluded for touching the image border.
     *
     * Expected output:
     * tuple(meta, regions_zarr, empty_stats_json)
     */
    ch_regions_inputs = CELLPOSE_FINAL.out.labels
        .map { meta, labels -> tuple(meta.id, meta, labels) }
        .join(PYMIF_MIP.out.zarr.map { meta, mip -> tuple(meta.id, mip) })
        .map { id, meta, labels, mip -> tuple(meta, labels, mip) }

    PYMIF_ENDOMORPH_REGIONS(ch_regions_inputs)

    /*
     * Measures per-cell morphology and V-cadherin/myosin/actin band
     * intensities against the MIP and the bands above, and adds the
     * image-level empty-area statistics to the per-image summary.
     *
     * Expected output:
     * tuple(meta, cells_csv, summary_csv)
     */
    ch_endomorph_inputs = PYMIF_MIP.out.zarr
        .map { meta, mip -> tuple(meta.id, meta, mip) }
        .join(CELLPOSE_FINAL.out.labels.map { meta, labels -> tuple(meta.id, labels) })
        .join(PYMIF_ENDOMORPH_REGIONS.out.zarr.map { meta, regions, stats -> tuple(meta.id, regions, stats) })
        .map { id, meta, mip, labels, regions, stats -> tuple(meta, mip, labels, regions, stats) }

    PYMIF_ENDOMORPH_FEATURES(ch_endomorph_inputs)

    /*
     * Merges the per-dataset CSVs into one per-cell CSV and one per-image
     * summary CSV, written to params.merged_output_dir (the directory of the
     * launch .sh script when started through template_run.sh).
     */
    def merged_dir = (params.merged_output_dir ?: workflow.launchDir).toString()
    if (!new File(merged_dir).isAbsolute()) {
        throw new IllegalArgumentException("--merged_output_dir must be an absolute path: got '${merged_dir}'")
    }

    PYMIF_ENDOMORPH_MERGE(
        PYMIF_ENDOMORPH_FEATURES.out.csv.map { meta, cells, summary -> cells }.collect(),
        PYMIF_ENDOMORPH_FEATURES.out.csv.map { meta, cells, summary -> summary }.collect(),
        merged_dir
    )

}
