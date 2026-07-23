process ILASTIK_SUBREGION_PIXEL_CLASS {
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
    def t_arg     = meta.t_index != null ? "--t-index ${meta.t_index}" : (meta.t != null ? "--t-index ${meta.t}" : "")
    def z_idx_arg = meta.z_index != null ? "--z-index ${meta.z_index}" : (meta.z != null ? "--z-index ${meta.z}" : "")
    def z_min_arg = meta.z_min   != null ? "--z-min ${meta.z_min}"     : ""
    def z_max_arg = meta.z_max   != null ? "--z-max ${meta.z_max}"     : ""
    def halo_arg  = meta.halo    != null ? "--halo ${meta.halo}"       : ""

    """
    ilastik_subregion_pixel_class.py \
        --input-zarr "${input_zarr}" \
        --output-zarr "${master_output_zarr}" \
        --project "${project}" \
        --y-min ${meta.y_min} \
        --y-max ${meta.y_max} \
        --x-min ${meta.x_min} \
        --x-max ${meta.x_max} \
        ${t_arg} \
        ${z_idx_arg} \
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
