process HYPERSTACK {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/biocontainers/bftools:8.0.0--hdfd78af_0"

    input:
    tuple val(meta), path(slices)

    output:
    output:
    tuple val(meta), path("*.ome.tif"), emit: hyperstack
    path "versions.yml"               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_segmented_stack"

    """
    bfconvert \\
        ${args} \\
        -nogroup \\
        "sample_${meta.sample}_t%t_z%z_c%c.tif" \\
        "${prefix}.ome.tif"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats: \$(bfconvert -version | head -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_segmented_stack"
    """
    touch ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats: \$(bfconvert -version | head -n 1)
    END_VERSIONS
    """
}
