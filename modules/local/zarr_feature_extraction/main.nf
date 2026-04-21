process ZARR_FEATURE_EXTRACTION {
    tag "${meta.id}"

    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_store)
    val intensity_comp
    val label_comp

    output:
    tuple val(meta), path("${meta.id}.features.parquet"), emit: parquet

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zarr_feature_extraction.py \
        -z ${zarr_store} \
        -o ${meta.id}.features.parquet \
        --intensity "${intensity_comp}" \
        --labels "${label_comp}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.features.parquet"
    """
}
