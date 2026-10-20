#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PYMIF_CONVERSION                        } from '../../modules/local/pymif_conversion/main'
include { PYMIF_MIP                               } from '../../modules/local/pymif_mip/main'
include { PYMIF_FEATURE_EXTRACTION                } from '../../modules/local/pymif_feature_extraction/main'
include { PYMIF_FILTER                            } from '../../modules/local/pymif_filter/main'
include { PYMIF_CLEAN_CHANNEL                     } from '../../modules/local/pymif_clean_channel/main'
include { CELLPOSE_SAM_ZARR_MIF as CELLPOSE_SEED  } from '../../modules/local/cellpose_sam_zarr_mif/main'
include { CELLPOSE_SAM_ZARR_MIF as CELLPOSE_FINAL } from '../../modules/local/cellpose_sam_zarr_mif/main'
include { PYMIF_ENDOMORPH_REGIONS                 } from '../../modules/local/pymif_endomorph_regions/main'
include { PYMIF_ENDOMORPH_FEATURES                } from '../../modules/local/pymif_endomorph_features/main'
include { PYMIF_ENDOMORPH_MERGE                   } from '../../modules/local/pymif_endomorph_merge/main'

workflow {
    def repoRootDir = workflow.projectDir.getParent().getParent()

    def resolveIntermediate = { row, String column, String datasetId ->
        Utils.normalizeOptionalText(row[column])
            ? Utils.resolveOutput(row, column, null, datasetId)
            : null
    }

    ch_cellpose_model = Channel.value(
        params.cellpose_model_path
            ? file(params.cellpose_model_path, checkIfExists: true)
            : file("${repoRootDir}/assets/NO_MODEL_FILE")
    )

    def noImageSentinel = file("${repoRootDir}/assets/NO_FILE")

    def row_index = new java.util.concurrent.atomic.AtomicInteger(0)

    pymif_inputs_ch = Channel
        .fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def index = row_index.getAndIncrement()
            def input_dataset = file(row.input, checkIfExists: true)

            def dataset_id = Utils.normalizeOptionalText(row.dataset_id)
                ?: "${input_dataset.simpleName}_${index + 1}"

            Utils.requireColumnValue(row, 'seed_channels')
            def seed_halo = Utils.parseOptionalFloat(row.seed_halo, 'seed_halo')

            def meta = row + [
                id                       : dataset_id,
                zarr_output              : resolveIntermediate(row, 'output', dataset_id),
                mip_output               : resolveIntermediate(row, 'mip_output', dataset_id),
                seed_output              : resolveIntermediate(row, 'seed_output', dataset_id),
                features_output          : resolveIntermediate(row, 'features_output', dataset_id),
                filter_output            : resolveIntermediate(row, 'filter_output', dataset_id),
                filter_csv_output        : resolveIntermediate(row, 'filter_csv_output', dataset_id),
                clean_output             : resolveIntermediate(row, 'clean_output', dataset_id),
                labels_output            : resolveIntermediate(row, 'labels_output', dataset_id),
                endomorph_regions_output : resolveIntermediate(row, 'endomorph_regions_output', dataset_id),
                endomorph_cells_output   : resolveIntermediate(row, 'endomorph_cells_output', dataset_id),
                endomorph_summary_output : resolveIntermediate(row, 'endomorph_summary_output', dataset_id),

                seed_segmentation_channels: Utils.normalizeOptionalText(row.seed_channels),
                seed_halo                 : seed_halo != null ? seed_halo : 20.0,
                seed_diameter             : Utils.parseOptionalInt(row.seed_diameter, 'seed_diameter'),
                seed_niter                : Utils.parseOptionalInt(row.seed_niter, 'seed_niter'),
                seed_flow_threshold       : Utils.parseOptionalFloat(row.seed_flow_threshold, 'seed_flow_threshold'),
                seed_cellprob_threshold   : Utils.parseOptionalFloat(row.seed_cellprob_threshold, 'seed_cellprob_threshold'),
                seed_use_physical_units   : Utils.normalizeOptionalText(row.seed_use_physical_units)?.toBoolean() ?: false,
                seed_do_3D                : Utils.normalizeOptionalText(row.seed_do_3D)?.toBoolean() ?: false,
            ]

            tuple(meta, input_dataset)
        }

    PYMIF_CONVERSION(pymif_inputs_ch)

    def boolOrFalse = { value ->
        def text = Utils.normalizeOptionalText(value)
        text ? text.toBoolean() : false
    }

    def cropBounds = { m ->
        [
            x_min: Utils.parseOptionalFloat(m.x_min, 'x_min'),
            x_max: Utils.parseOptionalFloat(m.x_max, 'x_max'),
            y_min: Utils.parseOptionalFloat(m.y_min, 'y_min'),
            y_max: Utils.parseOptionalFloat(m.y_max, 'y_max'),
            z_min: Utils.parseOptionalFloat(m.z_min, 'z_min'),
            z_max: Utils.parseOptionalFloat(m.z_max, 'z_max'),
        ]
    }

    def seedMeta = { m ->
        [
            id                   : m.id,
            zarr_output          : m.seed_output,
            segmentation_channels: m.seed_segmentation_channels,
            halo                 : m.seed_halo,
            diameter             : m.seed_diameter,
            niter                : m.seed_niter,
            flow_threshold       : m.seed_flow_threshold,
            cellprob_threshold   : m.seed_cellprob_threshold,
            use_physical_units   : m.seed_use_physical_units,
            do_3D                : m.seed_do_3D,
            resolution_level     : 0,
        ] + cropBounds(m)
    }

    def featuresMeta = { m ->
        def level = Utils.parseOptionalInt(m.features_resolution_level, 'features_resolution_level')
        [
            id                                  : m.id,
            csv_output                          : m.features_output,
            features_intensity_resolution_level : level != null ? level : 0,
            features_labels_resolution_level    : level != null ? level : 0,
            features_channels                   : Utils.normalizeOptionalText(m.features_channels),
            features_extra_properties           : Utils.normalizeOptionalText(m.features_extra_properties),
            features_extra_intensity_properties : Utils.normalizeOptionalText(m.features_extra_intensity_properties),
        ]
    }

    def filterMeta = { m ->
        [
            id              : m.id,
            zarr_output     : m.filter_output,
            csv_output      : m.filter_csv_output,
            filter_min_area : Utils.parseOptionalFloat(m.filter_min_area, 'filter_min_area'),
            filter_max_area : Utils.parseOptionalFloat(m.filter_max_area, 'filter_max_area'),
            filter_properties: Utils.normalizeOptionalText(m.filter_properties),
        ]
    }

    def cleanMeta = { m ->
        def clean_channel = Utils.parseOptionalInt(m.clean_channel, 'clean_channel')
        def default_clean_channel = (Utils.normalizeOptionalText(m.seed_channels).split(/[;,\s]+/).findAll { it })[0] as int
        [
            id           : m.id,
            zarr_output  : m.clean_output,
            clean_channel: clean_channel != null ? clean_channel : default_clean_channel,
        ]
    }

    def finalMeta = { m ->
        def level = Utils.parseOptionalInt(m.resolution_level, 'resolution_level')
        def halo  = Utils.parseOptionalFloat(m.cellpose_halo, 'cellpose_halo')
        [
            id                   : m.id,
            zarr_output          : m.labels_output,
            resolution_level     : level != null ? level : (params.resolution_level != null ? params.resolution_level as int : 0),
            use_physical_units   : boolOrFalse(m.use_physical_units),
            do_3D                : boolOrFalse(m.do_3D),
            diameter             : Utils.parseOptionalInt(m.cellpose_diameter, 'cellpose_diameter'),
            niter                : Utils.parseOptionalInt(m.cellpose_niter, 'cellpose_niter'),
            segmentation_channels: Utils.normalizeOptionalText(m.cellpose_channels),
            flow_threshold       : Utils.parseOptionalFloat(m.cellpose_flow_threshold, 'cellpose_flow_threshold'),
            cellprob_threshold   : Utils.parseOptionalFloat(m.cellpose_cellprob_threshold, 'cellpose_cellprob_threshold'),
            halo                 : halo != null ? halo : (params.cellpose_halo != null ? params.cellpose_halo as float : 20.0),
        ] + cropBounds(m)
    }

    def endomorphIntensityLevel = { m ->
        def level = Utils.parseOptionalInt(m.resolution_level, 'resolution_level')
        level != null ? level : (params.resolution_level != null ? params.resolution_level as int : 0)
    }

    def endomorphRingDepths = { m ->
        Utils.normalizeOptionalText(m.endomorph_ring_depths) ?: '6 9 24'
    }

    def endomorphTimepoint = { m ->
        def timepoint = Utils.parseOptionalInt(m.endomorph_timepoint, 'endomorph_timepoint')
        timepoint != null ? timepoint : 0
    }

    def floatOr = { m, String column, float fallback ->
        def value = Utils.parseOptionalFloat(m[column], column)
        value != null ? value : fallback
    }

    def regionsMeta = { m ->
        def border_margin = Utils.parseOptionalInt(m.endomorph_border_margin, 'endomorph_border_margin')
        [
            id                                  : m.id,
            zarr_output                         : m.endomorph_regions_output,
            endomorph_ring_depths               : endomorphRingDepths(m),
            endomorph_border_margin             : border_margin != null ? border_margin : 10,
            endomorph_timepoint                 : endomorphTimepoint(m),
            endomorph_intensity_resolution_level: endomorphIntensityLevel(m),
            endomorph_labels_resolution_level   : 0,
            endomorph_empty_channels            : Utils.normalizeOptionalText(m.endomorph_empty_channels),
            endomorph_empty_smoothing           : floatOr(m, 'endomorph_empty_smoothing', 2.0f),
            endomorph_empty_min_area            : floatOr(m, 'endomorph_empty_min_area', 25.0f),
            endomorph_empty_seed_percentile     : floatOr(m, 'endomorph_empty_seed_percentile', 1.0f),
            endomorph_empty_min_contrast        : floatOr(m, 'endomorph_empty_min_contrast', 0.25f),
        ]
    }

    def endomorphFeaturesMeta = { m ->
        [
            id                                  : m.id,
            csv_output                          : m.endomorph_cells_output,
            summary_csv_output                  : m.endomorph_summary_output,
            endomorph_ring_depths               : endomorphRingDepths(m),
            endomorph_timepoint                 : endomorphTimepoint(m),
            endomorph_intensity_resolution_level: endomorphIntensityLevel(m),
            endomorph_labels_resolution_level   : 0,
            endomorph_measure_channels          : Utils.normalizeOptionalText(m.endomorph_measure_channels),
        ]
    }

    PYMIF_MIP(
        PYMIF_CONVERSION.out.zarr.map { meta, zarr_path ->
            tuple(meta + [zarr_output: meta.mip_output], zarr_path)
        }
    )

    ch_mip = PYMIF_MIP.out.zarr

    ch_seed_inputs = ch_mip
        .map { meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(zarr_path, noImageSentinel)
            tuple(seedMeta(meta), local_image, image_uri)
        }

    CELLPOSE_SEED(ch_seed_inputs, ch_cellpose_model)

    ch_seed_labels = CELLPOSE_SEED.out.labels.map { meta, labels -> tuple(meta.id, labels) }

    ch_features_inputs = ch_mip
        .map { meta, zarr -> tuple(meta.id, featuresMeta(meta), zarr) }
        .join(ch_seed_labels)
        .map { id, features_meta, intensity, labels -> tuple(features_meta, intensity, labels) }

    PYMIF_FEATURE_EXTRACTION(ch_features_inputs)

    ch_filter_inputs = ch_mip
        .map { meta, zarr -> tuple(meta.id, filterMeta(meta)) }
        .join(ch_seed_labels)
        .join(PYMIF_FEATURE_EXTRACTION.out.csv.map { meta, csv -> tuple(meta.id, csv) })
        .map { id, filter_meta, labels, csv -> tuple(filter_meta, labels, csv) }

    PYMIF_FILTER(ch_filter_inputs)

    ch_clean_inputs = ch_mip
        .map { meta, intensity -> tuple(meta.id, cleanMeta(meta), intensity) }
        .join(PYMIF_FILTER.out.zarr_csv.map { meta, labels, csv -> tuple(meta.id, labels) })
        .map { id, clean_meta, intensity, labels -> tuple(clean_meta, intensity, labels) }

    PYMIF_CLEAN_CHANNEL(ch_clean_inputs)

    ch_final_inputs = ch_mip
        .map { meta, zarr -> tuple(meta.id, meta) }
        .join(PYMIF_CLEAN_CHANNEL.out.zarr.map { meta, zarr_path -> tuple(meta.id, zarr_path) })
        .map { id, meta, zarr_path ->
            def (local_image, image_uri) = Utils.splitImageInput(zarr_path, noImageSentinel)
            tuple(finalMeta(meta), local_image, image_uri)
        }

    CELLPOSE_FINAL(ch_final_inputs, ch_cellpose_model)

    ch_final_labels = CELLPOSE_FINAL.out.labels.map { meta, labels -> tuple(meta.id, labels) }

    ch_regions_inputs = ch_mip
        .map { meta, mip -> tuple(meta.id, regionsMeta(meta), mip) }
        .join(ch_final_labels)
        .map { id, regions_meta, mip, labels -> tuple(regions_meta, labels, mip) }

    PYMIF_ENDOMORPH_REGIONS(ch_regions_inputs)

    ch_endomorph_inputs = ch_mip
        .map { meta, mip -> tuple(meta.id, endomorphFeaturesMeta(meta), mip) }
        .join(ch_final_labels)
        .join(PYMIF_ENDOMORPH_REGIONS.out.zarr.map { meta, regions -> tuple(meta.id, regions) })
        .map { id, features_meta, mip, labels, regions -> tuple(features_meta, mip, labels, regions) }

    PYMIF_ENDOMORPH_FEATURES(ch_endomorph_inputs)

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
