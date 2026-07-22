process PYMIF_CONVERSION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.5.1'

    input:
    tuple val(meta), val(row), path(ome_tiff)

    output:
    tuple val(meta), path("*.zarr"), emit: zarr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // to prevent data duplication, 
    def local_row = new LinkedHashMap(row)
    local_row['input'] = ome_tiff.name

    // write to csv so batch2zarr can process
    def headers = local_row.keySet().join(',')
    def values  = local_row.values().join(',')
    def csv_content = """
    ${headers}
    ${values}
    """.stripIndent().trim()

    """
    cat <<'EOF' > input.csv
${csv_content}
EOF

    pymif batch2zarr \
        -i input.csv \
        ${args}

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
