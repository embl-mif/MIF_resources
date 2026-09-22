process PREALLOCATE_ZARR {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(input_zarr), path(tiles_csv), val(output_zarr)

    output:
    tuple val(meta), val(output_zarr), emit: zarr
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def num_channels_arg   = meta.num_channels   != null ? "--num-channels ${meta.num_channels}"                   : ""
    def channel_labels_arg = meta.channel_labels != null ? "--channel-labels \"${meta.channel_labels.join(',')}\"" : ""
    def dtype_arg          = meta.dtype          != null ? "--dtype ${meta.dtype}"                                 : ""
    def zarr_format_arg    = meta.zarr_format    != null ? "--zarr-format ${meta.zarr_format}"                     : ""
    def data_type_arg      = meta.data_type      != null ? "--data-type ${meta.data_type}"                         : ""
    def num_levels_arg     = meta.num_levels     != null ? "--num-levels ${meta.num_levels}"                       : ""

    // Crop bounds -- only meaningful alongside --input-zarr; preallocate_zarr.py
    // itself rejects these if passed without it. z/y/x also shift that axis'
    // translation; t only shrinks shape + records t_range (see its docstring).
    def t_min_arg = meta.t_min != null ? "--t-min ${meta.t_min}" : ""
    def t_max_arg = meta.t_max != null ? "--t-max ${meta.t_max}" : ""
    def z_min_arg = meta.z_min != null ? "--z-min ${meta.z_min}" : ""
    def z_max_arg = meta.z_max != null ? "--z-max ${meta.z_max}" : ""
    def y_min_arg = meta.y_min != null ? "--y-min ${meta.y_min}" : ""
    def y_max_arg = meta.y_max != null ? "--y-max ${meta.y_max}" : ""
    def x_min_arg = meta.x_min != null ? "--x-min ${meta.x_min}" : ""
    def x_max_arg = meta.x_max != null ? "--x-max ${meta.x_max}" : ""

    def from_tiles = tiles_csv.name != 'NO_FILE'

    if (!input_zarr) {
    error "PREALLOCATE_ZARR (${meta.id}): input_zarr (the original OME-Zarr image) is required."
}

    def command
    if (from_tiles) {
        command = """
        params=\$(calculate_zarr_params.py --tiles-csv "${tiles_csv}")
        preallocate_zarr.py \\
            --input-zarr "${input_zarr}" \\
            --output-zarr "${output_zarr}" \\
            --drop-axes c \\
            ${dtype_arg} \\
            ${zarr_format_arg} \\
            ${data_type_arg} \\
            ${num_levels_arg} \\
            \$params \\
            ${args}
        """
    } else {
        command = """
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
            ${z_min_arg} \\
            ${z_max_arg} \\
            ${y_min_arg} \\
            ${y_max_arg} \\
            ${x_min_arg} \\
            ${x_max_arg} \\
            ${num_levels_arg} \\
            ${args}
        """
    }

    """
    ${command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
        ngio: \$(python -c "import ngio; print(getattr(ngio, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
