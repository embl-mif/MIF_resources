process PYMIF_CONVERSION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.7.2'

    input:
    tuple val(meta), val(row), path(input_dataset)

    output:
    tuple val(meta), path("*.zarr"), emit: zarr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def isPresent = { value ->
        value != null &&
        value.toString().trim() &&
        !['-1', 'none', 'null', 'nan'].contains(value.toString().trim().toLowerCase())
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
        '--input_path', shellQuote(row.input),
        '--zarr_path', shellQuote(row.output_name)
    ]

    if (isPresent(row.microscope)) {
        commandParts.addAll(['--microscope', shellQuote(row.microscope)])
    }

    // CSV header "max_size(MB)" maps to --max_size. An explicit chunk size
    // takes precedence, matching PyMIF's batch conversion behavior.
    if (isPresent(row.chunk_size)) {
        commandParts.add('--chunk_size')
        commandParts.addAll(numericValues(row.chunk_size).collect { shellQuote(it) })
    } else if (isPresent(row['max_size(MB)'])) {
        commandParts.addAll(['--max_size', shellQuote(row['max_size(MB)'])])
    }

    [
        scene_index: '--scene_index',
        zarr_format: '--zarr_format',
        num_levels: '--num_levels',
        subset: '--subset'
    ].each { header, option ->
        if (isPresent(row[header])) {
            commandParts.addAll([option, shellQuote(row[header])])
        }
    }

    [
        downscale_factor: '--downscale_factor',
        channel_names: '--channel_names',
        channel_colors: '--channel_colors'
    ].each { header, option ->
        if (isPresent(row[header])) {
            def values = header == 'downscale_factor' ?
                numericValues(row[header]) :
                channelValues(row[header])
            commandParts.add(option)
            commandParts.addAll(values.collect { shellQuote(it) })
        }
    }

    def command = commandParts.join(' ')

    """
    ${command} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}.zarr"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: "stub"
    END_VERSIONS
    """
}
