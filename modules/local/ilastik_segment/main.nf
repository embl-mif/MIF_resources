process ILASTIK_SEGMENT {
    tag "${meta.id}"
    label 'process_medium'

    // To do: find container, using EMBL cluster module for now
    beforeScript 'module load ilastik/1.3.3post3'

    input:
    tuple val(meta), path(image)
    path  project

    output:
    tuple val(meta), path("*.tif"), emit: segmentation
    path "versions.yml"           , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}_segmented"
    def memory_mb = task.memory ? task.memory.toMega() : 4000

    """
    export LAZYFLOW_THREADS=${task.cpus}
    export LAZYFLOW_TOTAL_RAM_MB=${memory_mb}
    export HOME=\$PWD

    run_ilastik.sh \\
        --headless \\
        --project=${project} \\
        --raw_data ${image} \\
        ${args}

    mv *Simple*.tif ${prefix}.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(run_ilastik.sh --version 2>&1 | head -n 1 | sed 's/ilastik //')
    END_VERSIONS
    """
}
