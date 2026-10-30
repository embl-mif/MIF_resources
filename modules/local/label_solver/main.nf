process LABEL_SOLVER {
    tag "$meta.id"
    label "process_low"

    container "ghcr.io/embl-mif/multicut_stitching:0.8.0"

    input:
    tuple val(meta), path(local_graphs)

    output:
    tuple val(meta), path("*.json"), emit: global_mapping

    script:
    def beta_arg    = Utils.optionalCliArg("--beta", meta.beta)
    def output_json = "${meta.id}_global_label_mapping.json"

    """
    label_solver.py \\
        --local-graphs ${local_graphs} \\
        ${beta_arg} \\
        --output ${output_json}
    """

    stub:
    def output_json = "${meta.id}_global_label_mapping.json"
    """
    cat > ${output_json} <<-END_GLOBAL_MAPPING
    {}
    END_GLOBAL_MAPPING
    """
}

