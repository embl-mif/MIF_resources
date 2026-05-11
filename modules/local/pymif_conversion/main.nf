process PYMIF_CONVERSION {
    tag "${meta.id}"

    label 'process_medium'

    container 'ghcr.io/grinic/pymif:2026.5.1'

    input:
    tuple val(meta), path(input_csv)

    output:
    tuple val(meta), path("*.ome.zarr"), emit: zarr

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nextflowManagedArg = (task.ext.nextflow_managed != null ? (task.ext.nextflow_managed ? '--nextflow-managed' : '') : '--nextflow-managed')

    """
    echo "DEBUG: Running pymif batch2zarr for ID: ${meta.id}"
    pymif batch2zarr \
        -i ${input_csv} \
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.ome.zarr"
    """
}
