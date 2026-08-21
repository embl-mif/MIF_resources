process PREALLOCATE_ZARR {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(input_zarr), val(output_zarr)

    output:
    tuple val(meta), val(output_zarr), emit: empty_zarr
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def num_channels_arg   = meta.num_channels   != null ? "--num-channels ${meta.num_channels}"                     : ""
    def channel_labels_arg = meta.channel_labels != null ? "--channel-labels \"${meta.channel_labels.join(',')}\""   : ""
    def dtype_arg          = meta.dtype          != null ? "--dtype ${meta.dtype}"                                   : ""
    def zarr_format_arg    = meta.zarr_format    != null ? "--zarr-format ${meta.zarr_format}"                      : ""
    def data_type_arg      = meta.data_type      != null ? "--data-type ${meta.data_type}" : ""
    def t_min_arg          = meta.t_min          != null ? "--t-min ${meta.t_min}" : ""
    def t_max_arg          = meta.t_max          != null ? "--t-max ${meta.t_max}" : ""
    """
    preallocate_zarr.py \\
        --input-zarr "${input_zarr}" \\
        --output-zarr "${output_zarr}" \\
        ${num_channels_arg} \\
        ${channel_labels_arg} \\
        ${dtype_arg} \\
        ${zarr_format_arg} \\
        ${data_type_arg} \\
        ${t_min_arg} \\
        ${t_max_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
        ngio: \$(python -c "import ngio; print(getattr(ngio, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
