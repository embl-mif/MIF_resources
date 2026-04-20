process FEATURE_EXTRACTION {
    tag "${meta.id}"

    label 'process_medium'
    container "ghcr.io/grinic/pymif:2026.04.09"

    input:
    tuple val(meta), path(zarr, stageAs: 'input.zarr')

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    feature_extract.py \
        $args \
        -z input.zarr \
        -o "${prefix}.csv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.csv"
    """
}
