process FIND_TILE_NEIGHBORS {
    tag "$meta.id"
    label "process_low"

    container "docker://registry.git.embl.org/grp-cba/containers/zarr_tools:0.34.0"

    input:
    tuple val(meta), path(tiles_csv)

    output:
    tuple val(meta), path("*.csv"), emit: tile_neighbors

    script:
    def output_csv = "${meta.id}_tile_neighbors.csv"
    """
    tile_neighbors.py \\
        --tiles-csv ${tiles_csv} \\
        --output ${output_csv}
    """

    stub:
    def output_csv = "${meta.id}_tile_neighbors.csv"
    """
    cat > ${output_csv} <<-END_NEIGHBORS_CSV
    source_dataset_id,dataset_id,segmented_tile_path,right_dataset_id,right_tile_index,right_segmented_tile_path,bottom_dataset_id,bottom_tile_index,bottom_segmented_tile_path,behind_dataset_id,behind_tile_index,behind_segmented_tile_path
    ${meta.id},${meta.id}__tile_000001,tile_000001_labels.ome.zarr,,,,,,,,,
    END_NEIGHBORS_CSV
    """
}
