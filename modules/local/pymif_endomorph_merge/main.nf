process PYMIF_ENDOMORPH_MERGE {
    label 'process_low'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    path cells_csvs, stageAs: 'cells/*'
    path summary_csvs, stageAs: 'summaries/*'
    val output_cells
    val output_summary

    output:
    path "versions.yml", emit: versions

    script:
    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    """
    mkdir -p \$(dirname ${shellQuote(output_cells)})
    mkdir -p \$(dirname ${shellQuote(output_summary)})

    merge_endomorph_csvs.py \
        --cells cells/* \
        --summaries summaries/* \
        --output-cells ${shellQuote(output_cells)} \
        --output-summary ${shellQuote(output_summary)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
