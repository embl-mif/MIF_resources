process WRITE_RELABELED_TILE {
    tag "$meta.id"
    label "process_low"

    container "docker://registry.git.embl.org/grp-cba/containers/zarr_tools:0.34.0"

    input:
    tuple val(meta), path(tile_label), path(tiles_csv), path(store), path(global_mapping)

    output:
    tuple val(meta), path("*.done"), emit: write_confirmation

    script:
    def output_marker = "${meta.id}.done"

    if (meta.label_offset_step == null) {
        throw new IllegalArgumentException("label_offset_step is required in meta for WRITE_RELABELED_TILE")
    }

    def z_min_arg         = Utils.optionalCliArg("--z-min", meta.z_min)
    def tile_size_z_arg   = Utils.optionalCliArg("--tile-size-z", meta.tile_size_z)
    def tile_overlap_z_arg = Utils.optionalCliArg("--tile-overlap-z", meta.tile_overlap_z)

    // Whether this tile has a real forward neighbor on each axis -- decides whether its own
    // core region (the part it writes into the shared store) gets trimmed on that axis, see
    // bin/write_relabeled_tile.py. Same right/bottom/behind_dataset_id fields FIND_TILE_NEIGHBORS
    // already produces and LOCAL_LABEL_MATCHER already consumes for a different purpose.
    def has_right_arg  = meta.right_dataset_id  ? "--has-right"  : ""
    def has_bottom_arg = meta.bottom_dataset_id ? "--has-bottom" : ""
    def has_behind_arg = meta.behind_dataset_id ? "--has-behind" : ""

    """
    write_relabeled_tile.py \\
        --tile-label-path ${tile_label} \\
        --tiles-csv ${tiles_csv} \\
        --store-path ${store} \\
        --global-mapping ${global_mapping} \\
        --tile-index ${meta.tile_index} \\
        --label-offset-step ${meta.label_offset_step} \\
        --x-min ${meta.x_min} \\
        --y-min ${meta.y_min} \\
        ${z_min_arg} \\
        --tile-size-x ${meta.tile_size_x} \\
        --tile-size-y ${meta.tile_size_y} \\
        ${tile_size_z_arg} \\
        --tile-overlap-x ${meta.tile_overlap_x} \\
        --tile-overlap-y ${meta.tile_overlap_y} \\
        ${tile_overlap_z_arg} \\
        ${has_right_arg} \\
        ${has_bottom_arg} \\
        ${has_behind_arg} \\
        --output ${output_marker}
    """

    stub:
    def output_marker = "${meta.id}.done"
    """
    touch ${output_marker}
    """
}

