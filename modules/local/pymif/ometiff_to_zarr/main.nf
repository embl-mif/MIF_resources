process CONVERT_TO_ZARR {
    tag "${meta.exp_name}-${meta.well}-${meta.field}"
    label 'process_medium'
    container "ghcr.io/grinic/pymif:2026.04.09"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.zarr"), emit: zarr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.well}_${meta.field}"

    def nl_arg = meta.num_layers ? "-nl ${meta.num_layers}" : ""
    def cs_arg = meta.chunk_size ? "-cs '${meta.chunk_size.join(' ')}'" : ""

    """
    ometiff2zarr.py \\
        -i $tiff \\
        -o ${prefix}.zarr \\
        $nl_arg \\
        $cs_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
