process ILASTIK_ZARR {
    tag "${meta.id}_z${meta.z_min}_y${meta.y_min}"
    label 'process_medium'

    container "docker.io/biocontainers/ilastik:1.4.2_cv1"

    input:
    tuple val(meta), val(input_zarr), val(master_output_zarr)
    path  project

    output:
    tuple val(meta), val(master_output_zarr), emit: zarr
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    // Dynamically handle time, 2D slice indices, 3D bounds, and halo from meta map
    def t_min_arg = meta.t_min != null ? "--t-min ${meta.t_min}" : ""
    def t_max_arg = meta.t_max != null ? "--t-max ${meta.t_max}" : ""
    def y_min_arg = meta.y_min   != null ? "--y-min ${meta.y_min}"     : ""
    def y_max_arg = meta.y_max   != null ? "--y-max ${meta.y_max}"     : ""
    def x_min_arg = meta.x_min   != null ? "--x-min ${meta.x_min}"     : ""
    def x_max_arg = meta.x_max   != null ? "--x-max ${meta.x_max}"     : ""
    def z_min_arg = meta.z_min   != null ? "--z-min ${meta.z_min}"     : ""
    def z_max_arg = meta.z_max   != null ? "--z-max ${meta.z_max}"     : ""
    def halo_arg  = meta.halo    != null ? "--halo ${meta.halo}"       : ""

    def ilastik_python = "/opt/ilastik-1.4.2-Linux/bin/python3"

    """
    LAZYFLOW_THREADS=${task.cpus} LAZYFLOW_TOTAL_RAM_MB=${task.memory.toMega()} \\
    ${ilastik_python} ${moduleDir}/resources/usr/bin/ilastik_tile.py \
        --input-zarr "${input_zarr}" \
        --output-zarr "${master_output_zarr}" \
        --project "${project}" \
        ${y_min_arg} \
        ${y_max_arg} \
        ${x_min_arg} \
        ${x_max_arg} \
        ${t_min_arg} \
        ${t_max_arg} \
        ${z_min_arg} \
        ${z_max_arg} \
        ${halo_arg} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(python -c "import ilastik; print(getattr(ilastik, '__version__', '1.4.2'))")
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
        h5py: \$(python -c "import h5py; print(h5py.__version__)")
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: 1.4.2
    END_VERSIONS
    """
}
