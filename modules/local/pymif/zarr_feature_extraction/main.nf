process FEATURE_EXTRACTION {
    tag "${meta.exp_name}-${meta.well}-${meta.field}"

    label 'process_medium'
    container "ghcr.io/grinic/pymif:2026.04.09"

    input:
    tuple val(meta), path(zarr, stageAs: 'input.zarr')

    output:
    tuple val(meta), path("*.csv"), emit: csv

    script:
    """
    feature_extract.py \
        -z input.zarr \
        -o "${meta.well}_${meta.field}.csv"
    """
}
