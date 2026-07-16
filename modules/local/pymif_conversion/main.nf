process PYMIF_CONVERSION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.5.1'

    input:
    tuple val(meta), path(row)

    output:
    tuple val(meta), path("*.zarr"), emit: zarr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def software = 'pymif'

    // write to csv so batch2zarr can process
    def headers = row.keySet().join(',')
    def values  = row."""
    ${headers}
    ${values}
    """.stripIndent().trim()

    """
    cat <<'EOF' > input.csv
${csv_content}
EOF

    pymif batch2zarr \
        -i input.csv \
        -o ${prefix}.zarr \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${software}: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

}
