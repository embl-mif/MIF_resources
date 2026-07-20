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

    // Handle optional Z flags dynamically for 2D vs 3D data
    def z_min_arg = meta.z_min != null ? "--z-min ${meta.z_min}" : ""
    def z_max_arg = meta.z_max != null ? "--z-max ${meta.z_max}" : ""
    def halo_arg  = meta.halo  != null ? "--halo ${meta.halo}"   : ""

    """
    ilastik_subregion_pixel_class.py \
        --input-zarr "${input_zarr}" \
        --output-zarr "${master_output_zarr}" \
        --project "${project}" \
        --y-min ${meta.y_min} \
        --y-max ${meta.y_max} \
        --x-min ${meta.x_min} \
        --x-max ${meta.x_max} \
        ${z_min_arg} \
        ${z_max_arg} \
        ${halo_arg} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ilastik: \$(python -c "import ilastik; print(getattr(ilastik, '__version__', '1.4.2'))")
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
