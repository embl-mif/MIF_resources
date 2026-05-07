process ZARR_FEATURE_EXTRACTION {
    tag "${meta.id}"

    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(zarr_store)
    val intensity_comp
    val label_comp

    output:
    tuple val(meta), path("*.parquet"), emit: parquet

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nextflowManagedArg = (task.ext.nextflow_managed != null ? (task.ext.nextflow_managed ? '--nextflow-managed' : '') : '--nextflow-managed')

    """
    export MODULE_LIB="${moduleDir}/resources/usr/bin/lib"
    extract_features.py \
        -z ${zarr_store} \
        -o ${prefix}.features.parquet \
        --intensity "${intensity_comp}" \
        --segmentation "${label_comp}" \
        ${nextflowManagedArg} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.features.parquet"
    """
}
