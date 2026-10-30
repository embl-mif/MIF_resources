process PREALLOCATE_ZARR {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(local_input_zarr), val(input_zarr_uri), path(tiles_csv)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}.ome.zarr"
    def input_zarr = input_zarr_uri ?: local_input_zarr

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

    def from_tiles = !tiles_csv.name.startsWith('NO_')

    // local_input_zarr is always a real, non-null staged path -- either the actual local
    // store or a NO_*-prefixed sentinel -- so presence must be checked by the sentinel
    // naming convention, not truthiness or an exact 'NO_FILE' match, since different
    // callers use differently-named sentinels for this slot (NO_FILE, NO_RAW_IMAGE, etc.)
    // to avoid Nextflow input-file-collision errors when multiple path slots in one call
    // are simultaneously absent.
    def hasLocalInput = !local_input_zarr.name.startsWith('NO_')
    if (!input_zarr_uri && !hasLocalInput) {
        error "PREALLOCATE_ZARR (${meta.id}): input_zarr (the original OME-Zarr image) is required."
    }

    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    // NOTE: no `def` here -- a `def`-scoped variable is invisible to the output: block's
    // own evaluation, so path(output_zarr_link) would look for a literal file named
    // "output_zarr_link" instead of resolving to this value (see ZARRSMITH_CONVERSION for
    // the same pattern done correctly).
    output_zarr_link = new File(output_zarr.toString()).name

    def command
    if (from_tiles) {
        command = """
        params=\$(calculate_zarr_params.py --tiles-csv "${tiles_csv}" --source-zarr "${input_zarr}")
        preallocate_zarr.py \\
            --input-zarr "${input_zarr}" \\
            --output-zarr "${output_zarr}" \\
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

    ${Utils.shellLinkOutput(output_zarr)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
        ngio: \$(python -c "import ngio; print(getattr(ngio, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
