process SLICE_IMAGE {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/biocontainers/bftools:8.0.0--hdfd78af_0"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.tif"), emit: single_slices
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    // Default naming pattern if not overridden in config
    def prefix = task.ext.prefix ?: "${meta.id}_slice"
    
    """
    # bfconvert uses %t for timepoint, %z for z-slice, and %c for channel indices
    bfconvert \\
        ${args} \\
        "${image}" \\
        "${prefix}_t%t_z%z_c%c.tif"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioformats: \$(bfconvert -version | head -n 1)
    END_VERSIONS
    """
}
