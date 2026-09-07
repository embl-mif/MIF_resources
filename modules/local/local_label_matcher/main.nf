process LOCAL_LABEL_MATCHER {
    tag "$meta.id"
    label "process low"

    container "docker://registry.git.embl.org/grp-cba/containers/multicut_stitching:0.8.0"

    input:
    tuple val(meta), path(self_zarr), path(right_zarr, stageAs: 'right_input'), path(bottom_zarr, stageAs: 'bottom_input'), path(behind_zarr, stageAs: 'behind_input')

    output:
    tuple val(meta), path("*.json"), emit: local_graph

    script:
    // right_zarr/bottom_zarr/behind_zarr are the NO_FILE sentinel (see assets/NO_FILE) when
    // that neighbor doesn't exist -- staged like any other input so this and every other
    // zarr path here is a real Nextflow `path`, auto-bind-mounted under container executors
    // instead of relying on an implicit shared filesystem (workflow_status.md). Presence is
    // gated on meta.*_dataset_id (blank when absent), not on inspecting the staged path
    // itself -- `stageAs` above renames every input to its own fixed local name, so the
    // sentinel's *staged* name is indistinguishable from a real neighbor's.
    def right_path_arg   = meta.right_dataset_id  ? "--right-path ${right_zarr}"   : ""
    def right_index_arg  = Utils.optionalCliArg("--right-tile-index", meta.right_tile_index)
    def bottom_path_arg  = meta.bottom_dataset_id ? "--bottom-path ${bottom_zarr}" : ""
    def bottom_index_arg = Utils.optionalCliArg("--bottom-tile-index", meta.bottom_tile_index)
    def behind_path_arg  = meta.behind_dataset_id ? "--behind-path ${behind_zarr}" : ""
    def behind_index_arg = Utils.optionalCliArg("--behind-tile-index", meta.behind_tile_index)
    def z_overlap_arg    = Utils.optionalCliArg("--tile-overlap-z", meta.tile_overlap_z)
    def output_json      = "${meta.id}_local_graph.json"

    if (meta.label_offset_step == null) {
        throw new IllegalArgumentException("label_offset_step is required in meta for LOCAL_LABEL_MATCHER")
    }

    """
    label_matching.py \\
        --self-path ${self_zarr} \\
        --tile-index ${meta.tile_index} \\
        --label-offset-step ${meta.label_offset_step} \\
        --tile-overlap-x ${meta.tile_overlap_x} \\
        --tile-overlap-y ${meta.tile_overlap_y} \\
        ${z_overlap_arg} \\
        ${right_path_arg} ${right_index_arg} \\
        ${bottom_path_arg} ${bottom_index_arg} \\
        ${behind_path_arg} ${behind_index_arg} \\
        --output ${output_json}
    """

    stub:
    def output_json = "${meta.id}_local_graph.json"
    """
    cat > ${output_json} <<-END_LOCAL_GRAPH
    {"edges": [], "disaffinities": []}
    END_LOCAL_GRAPH
    """
}

