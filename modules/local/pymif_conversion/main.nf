process PYMIF_CONVERSION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.7.4'

    input:
    tuple val(meta), path(input_dataset)

    output:
    tuple val(meta), val(meta.zarr_output_name), emit: zarr
    path "versions.yml"                  , emit: versions

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
        '--input_path', shellQuote(meta.input),
        '--zarr_path', shellQuote(meta.zarr_output_name)
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
        subset: '--subset'
    ].each { header, option ->
        if (isPresent(meta[header])) {
            commandParts.addAll([option, shellQuote(meta[header])])
        }
    }

    [
        downscale_factor: '--downscale_factor',
        channel_names: '--channel_names',
        channel_colors: '--channel_colors'
    ].each { header, option ->
        if (isPresent(meta[header])) {
            def values = header == 'downscale_factor' ?
                numericValues(meta[header]) :
                channelValues(meta[header])
            commandParts.add(option)
            commandParts.addAll(values.collect { shellQuote(it) })
        }
    }

    def command = commandParts.join(' ')

    """
    if [ -d ${shellQuote(meta.zarr_output_name)} ] && [ -e ${shellQuote(meta.zarr_output_name)}/zarr.json -o -e ${shellQuote(meta.zarr_output_name)}/.zgroup ]; then
        echo "Output already exists, skipping conversion: ${meta.zarr_output_name}"
    else
        mkdir -p \$(dirname ${shellQuote(meta.zarr_output_name)})
        ${command} ${args}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.zarr_output_name}"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: "stub"
    END_VERSIONS
    """
}
