include { TILE_CONSTRUCTOR      } from '../../../modules/local/tileconstructor/main'
include { CELLPOSE_SAM_ZARR_MIF } from '../../../modules/local/cellpose_sam_zarr_mif/main'
include { PREALLOCATE_ZARR      } from '../../../modules/local/preallocate_zarr/main'
include { FIND_TILE_NEIGHBORS   } from '../../../modules/local/find_tile_neighbors/main'
include { LOCAL_LABEL_MATCHER   } from '../../../modules/local/local_label_matcher/main'
include { LABEL_SOLVER          } from '../../../modules/local/label_solver/main'
include { WRITE_RELABELED_TILE  } from '../../../modules/local/write_relabeled_tile/main'

workflow TILED_CELLPOSE_SAM_ZARR {

    take:
    ch_images              // channel: [ meta, image_path ] -- meta already carries
                            //   per-image segmentation_channels/halo/use_physical_units
    val_output_dir         // string:  base output dir; merged stores land under
                            //   "${val_output_dir}/<source_id>_merged_mask.ome.zarr"
    val_label_offset_step  // integer: per-tile global label-id offset step
    val_multicut_beta      // float:   LABEL_SOLVER multicut beta parameter

    main:

    ch_versions = Channel.empty()


    def noNeighborZarr  = file("${workflow.projectDir}/assets/NO_FILE")
    def noImageSentinel = file("${workflow.projectDir}/assets/NO_FILE")

    ch_split_images = ch_images.map { meta, image ->
        def (localImage, imageUri) = Utils.splitImageInput(image, noImageSentinel)
        tuple(meta.id, meta, image, localImage, imageUri)
    }

    TILE_CONSTRUCTOR( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(meta, imageUri ?: localImage.toString()) } )

    ch_merged_zarr_input = TILE_CONSTRUCTOR.out.tiles
        .map { parentMeta, tileCsv -> tuple(parentMeta.id, tileCsv) }
        .join( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(id, localImage, imageUri) } )
        .map { srcId, tileCsv, localImage, imageUri ->
            def sourceImagePath = imageUri ?: localImage.toString()
            def mergedOutputZarr = "${val_output_dir}/${srcId}_merged_mask.ome.zarr"
            def mergedMeta = [id: srcId, dtype: 'uint32', data_type: 'label']
            tuple(mergedMeta, '', tileCsv, sourceImagePath, mergedOutputZarr)
        }

    PREALLOCATE_ZARR(ch_merged_zarr_input)
    ch_versions = ch_versions.mix(PREALLOCATE_ZARR.out.versions)

    ch_image_tiles = TILE_CONSTRUCTOR.out.tiles
        .map { parentMeta, tileCsv -> tuple(parentMeta.id, tileCsv) }
        .join( ch_split_images.map { id, meta, image, localImage, imageUri -> tuple(id, meta, image, localImage, imageUri) } )
        .flatMap { _srcId, tileCsv, parentMeta, image, localImage, imageUri ->
            tileCsv.splitCsv(header: true).collect { row ->
                def tileMeta = [
                    id: row.dataset_id,
                    source_dataset_id: row.source_dataset_id,
                    tile_id: row.tile_id,
                    tile_index: row.tile_index,
                    resolution_level: row.resolution_level,
                    x_min: row.x_min,
                    x_max: row.x_max,
                    y_min: row.y_min,
                    y_max: row.y_max,
                    z_min: row.z_min,
                    z_max: row.z_max,
                    tile_overlap_x: row.tile_overlap_x,
                    tile_overlap_y: row.tile_overlap_y,
                    tile_overlap_z: row.tile_overlap_z,
                    tile_size_x: row.tile_size_x,
                    tile_size_y: row.tile_size_y,
                    tile_size_z: row.tile_size_z,
                    segmentation_channels: parentMeta.segmentation_channels,
                    halo: parentMeta.halo,
                    use_physical_units: false,
                    do_3D: parentMeta.do_3D ?: false,
                ]
 
                tuple(tileMeta, image, localImage, imageUri)
            }
        }

    ch_cellpose_model = params.cellpose_model_path
        ? Channel.value(file(params.cellpose_model_path, checkIfExists: true))
        : Channel.value(file("${workflow.projectDir}/assets/NO_MODEL_FILE"))
 
    CELLPOSE_SAM_ZARR_MIF(
        ch_image_tiles.map { meta, image, localImage, imageUri -> tuple(meta, localImage, imageUri) },
        ch_cellpose_model
    )

    ch_tiles_and_segmentations_csv = CELLPOSE_SAM_ZARR_MIF.out.labels
        .map { meta, label ->
            tuple(meta.source_dataset_id, "${meta.id},${meta.source_dataset_id},${meta.tile_id},${meta.tile_index},${meta.x_min},${meta.x_max},${meta.y_min},${meta.y_max},${meta.z_min},${meta.z_max},${meta.resolution_level},${meta.tile_overlap_x},${meta.tile_overlap_y},${meta.tile_overlap_z},${label}")
        }
        .collectFile(seed: 'dataset_id,source_dataset_id,tile_id,tile_index,x_min,x_max,y_min,y_max,z_min,z_max,resolution_level,tile_overlap_x,tile_overlap_y,tile_overlap_z,segmented_tile_path', newLine: true) { srcId, row ->
            [ "${srcId}_tiles_and_segmentations.csv", row ]
        }
        .map { csvPath -> tuple(csvPath.baseName.replace('_tiles_and_segmentations', ''), csvPath) }
 
    ch_tile_neighbors_input = ch_tiles_and_segmentations_csv
        .map { srcId, csv -> tuple([id: srcId], csv) }
 
    FIND_TILE_NEIGHBORS(ch_tile_neighbors_input)

    // dataset_id -> label path, for resolving neighbor labels below.
    ch_label_path_by_id = CELLPOSE_SAM_ZARR_MIF.out.labels
        .map { meta, label -> tuple(meta.id, label) }
        .toList()
        .map { pairs -> pairs.collectEntries { id, label -> [(id): label] } }
 
    ch_local_matcher_input = CELLPOSE_SAM_ZARR_MIF.out.labels
        .map { meta, label -> tuple(meta.id, meta + [label_offset_step: val_label_offset_step], label) }
        .join(
            FIND_TILE_NEIGHBORS.out.tile_neighbors.flatMap { srcMeta, csv ->
                csv.splitCsv(header: true).collect { row -> tuple(row.dataset_id, row) }
            }
        )
        .combine(ch_label_path_by_id)
        .map { tileId, meta, selfLabel, neighborRow, labelPathById ->
            def fullMeta = meta + [
                right_dataset_id : neighborRow.right_dataset_id,
                right_tile_index : neighborRow.right_tile_index,
                bottom_dataset_id: neighborRow.bottom_dataset_id,
                bottom_tile_index: neighborRow.bottom_tile_index,
                behind_dataset_id: neighborRow.behind_dataset_id,
                behind_tile_index: neighborRow.behind_tile_index,
            ]
 
            def rightZarr  = labelPathById[neighborRow.right_dataset_id]  ?: noNeighborZarr
            def bottomZarr = labelPathById[neighborRow.bottom_dataset_id] ?: noNeighborZarr
            def behindZarr = labelPathById[neighborRow.behind_dataset_id] ?: noNeighborZarr
 
            tuple(fullMeta, selfLabel, rightZarr, bottomZarr, behindZarr)
        }
 
    LOCAL_LABEL_MATCHER(ch_local_matcher_input)

 
    ch_local_graphs_grouped = LOCAL_LABEL_MATCHER.out.local_graph
        .map { meta, json -> tuple(meta.source_dataset_id, json) }
        .groupTuple()
        .map { srcId, jsonList -> tuple([id: srcId, beta: val_multicut_beta], jsonList) }
 
    LABEL_SOLVER(ch_local_graphs_grouped)
 
    ch_write_tile_input = CELLPOSE_SAM_ZARR_MIF.out.labels
        .map { meta, label -> tuple(meta.id, meta + [label_offset_step: val_label_offset_step], label) }
        .join(
            FIND_TILE_NEIGHBORS.out.tile_neighbors.flatMap { _srcMeta, csv ->
                csv.splitCsv(header: true).collect { row -> tuple(row.dataset_id, row) }
            }
        )
        .map { _tileId, meta, label, neighborRow ->
            def fullMeta = meta + [
                right_dataset_id : neighborRow.right_dataset_id,
                bottom_dataset_id: neighborRow.bottom_dataset_id,
                behind_dataset_id: neighborRow.behind_dataset_id,
            ]
            tuple(meta.source_dataset_id, fullMeta, label)
        }
        .combine( TILE_CONSTRUCTOR.out.tiles.map { parentMeta, tileCsv -> tuple(parentMeta.id, tileCsv) }, by: 0 )
        .combine( PREALLOCATE_ZARR.out.zarr.map { meta, store -> tuple(meta.id, store) }, by: 0 )
        .combine( LABEL_SOLVER.out.global_mapping.map { meta, json -> tuple(meta.id, json) }, by: 0 )
        .map { _srcId, fullMeta, label, tileCsv, store, mapping ->
            tuple(fullMeta, label, tileCsv, store, mapping)
        }
 
    WRITE_RELABELED_TILE(ch_write_tile_input)

 
    ch_write_confirmed = WRITE_RELABELED_TILE.out.write_confirmation
        .map { meta, confirmation -> tuple(meta.source_dataset_id, confirmation) }
        .groupTuple(by: 0)
 
    ch_mask_ready = PREALLOCATE_ZARR.out.zarr
        .map { meta, store -> tuple(meta.id, meta, store) }
        .join( ch_write_confirmed )
        .map { _srcId, meta, store, _confirmations -> tuple(meta, store) }
 
    emit:
    mask     = ch_mask_ready   // channel: [ meta, merged_mask.ome.zarr ]
    versions = ch_versions     // channel: [ versions.yml ]
}
