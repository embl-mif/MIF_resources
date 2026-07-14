process ILASTIK_SUBREGION_PIXEL_CLASS {
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/biocontainers/ilastik:1.4.2_cv1"

    input:
    tuple val(meta), val(image), val(output_zarr)
    path  project

    output:
    tuple val(meta), val(output_zarr), emit: segmentation
    path "versions.yml"           , emit: versions

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}_segmented"
    def memory_mb = task.memory ? task.memory.toMega() : 4000

    def subregion = ""
    def axes = "yxc"
    if (meta.y_min != null && meta.y_min != "") {
        if (meta.z_min != null && meta.z_min != "") {
            axes = "zyxc"
            subregion = "--cutout_subregion=\"[(${meta.z_min},${meta.y_min},${meta.x_min},0),(${meta.z_max},${meta.y_max},${meta.x_max},100)]\""
        } else {
            subregion = "--cutout_subregion=\"[(${meta.y_min},${meta.x_min},0),(${meta.y_max},${meta.x_max},100)]\""
        }
    }

    """
    export LAZYFLOW_THREADS=${task.cpus}
    export LAZYFLOW_TOTAL_RAM_MB=${memory_mb}
    export HOME=\$PWD

    run_ilastik.sh \\
        --headless \\
        --project=${project} \\
        --input_axes=${axes} \\
        --output_format="zarr" \\
        --output_filename_format="${output_zarr}" \\
        --output_internal_path="0" \\
        --export_source="probabilities" \\
        ${subregion} \\
        ${args} \\
        "${image}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(run_ilastik.sh --headless --version 2>&1 | head -n 1 | sed 's/ilastik //')
    END_VERSIONS
    """
}
