process PYMIF_CONVERSION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.10.0'

    input:
    tuple val(meta), path(input_dataset)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_zarr_link = new File(output_zarr.toString()).name

	def args = task.ext.args ?: ''

    def isPresent = { value ->
        value != null &&
        value.toString().trim() &&
        !['-1', 'none', 'null', 'nan'].contains(value.toString().trim().toLowerCase())
    }
    // "none" is a meaningful value for shard_exclude_axes (it disables the
    // default t/c exclusion so every axis may be merged into a shard), unlike
    // the other columns where "none" means "not set".
    def isPresentAllowNone = { value ->
        value != null &&
        value.toString().trim() &&
        !['-1', 'null', 'nan'].contains(value.toString().trim().toLowerCase())
    }
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }
    def numericValues = { value ->
        value.toString().trim().split(/[\s,;]+/).findAll { it }
    }
    def channelValues = { value ->
        def separator = value.toString().contains(',') ? /,\s*/ : /\s+/
        value.toString().trim().split(separator).findAll { it }
    }

    // Required CSV columns have CLI names that differ from their headers.
    def commandParts = [
        'pymif', '2zarr',
        '--input_path', shellQuote(input_dataset),
        '--zarr_path', shellQuote(output_zarr)
    ]

    if (isPresent(meta.microscope)) {
        commandParts.addAll(['--microscope', shellQuote(meta.microscope)])
    }

    // CSV header "max_size(MB)" maps to --max_size. An explicit chunk size
    // takes precedence, matching PyMIF's batch conversion behavior.
    if (isPresent(meta.chunk_size)) {
        commandParts.add('--chunk_size')
        commandParts.addAll(numericValues(meta.chunk_size).collect { shellQuote(it) })
    } else if (isPresent(meta['max_size(MB)'])) {
        commandParts.addAll(['--max_size', shellQuote(meta['max_size(MB)'])])
    }

    [
        scene_index: '--scene_index',
        zarr_format: '--zarr_format',
        num_levels: '--num_levels',
        subset: '--subset',
        shard_target_mb: '--shard_target_mb',
        num_workers: '--num_workers'
    ].each { header, option ->
        if (isPresent(meta[header])) {
            commandParts.addAll([option, shellQuote(meta[header])])
        }
    }

    [
        downscale_factor: '--downscale_factor',
        channel_names: '--channel_names',
        channel_colors: '--channel_colors',
        shards: '--shards'
    ].each { header, option ->
        // "none" is meaningful for shards (disables sharding; pymif defaults to auto).
        if (header == 'shards' ? isPresentAllowNone(meta[header]) : isPresent(meta[header])) {
            def values = header == 'downscale_factor' ?
                numericValues(meta[header]) :
                channelValues(meta[header])
            commandParts.add(option)
            commandParts.addAll(values.collect { shellQuote(it) })
        }
    }

    // shard_exclude_axes only applies with --shards, but is parsed on its own
    // so an explicit "none" (allow every axis to merge) is not dropped.
    if (isPresentAllowNone(meta.shard_exclude_axes)) {
        commandParts.add('--shard_exclude_axes')
        commandParts.addAll(channelValues(meta.shard_exclude_axes).collect { shellQuote(it) })
    }

    // pymif >= 2026.10 drops singleton t/c/z axes by default. Downstream modules expect
    // the full axis layout, so keep them unless the CSV explicitly sets drop_singleton.
    def dropSingleton = isPresent(meta.drop_singleton) ? meta.drop_singleton.toString().trim().toLowerCase() : 'false'
    commandParts.addAll(['--drop_singleton', shellQuote(dropSingleton)])

    def command = commandParts.join(' ')

    """
    if [ -d ${shellQuote(output_zarr)} ] && [ -e ${shellQuote(output_zarr)}/zarr.json -o -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "Output already exists, skipping conversion: ${output_zarr}"
    else
        mkdir -p \$(dirname ${shellQuote(output_zarr)})
        ${command} ${args}
    fi

    if [ ! -e ${shellQuote(output_zarr)}/zarr.json ] && [ ! -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_zarr}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_zarr)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    output_zarr = meta.zarr_output ?: "${meta.id}.ome.zarr"
    output_zarr_link = new File(output_zarr.toString()).name
    """
    mkdir -p "${output_zarr}"
    ${Utils.shellLinkOutput(output_zarr)}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: "stub"
    END_VERSIONS
    """
}
