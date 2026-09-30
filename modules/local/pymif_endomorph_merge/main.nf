process PYMIF_ENDOMORPH_MERGE {
    label 'process_low'

    publishDir { merged_dir }, mode: 'copy'

    input:
    path cells_csvs, stageAs: 'cells/?/*'
    path summary_csvs, stageAs: 'summaries/?/*'
    val merged_dir

    output:
    path "endomorph_cells_merged.csv"
    path "endomorph_summary_merged.csv"

    script:
    // Keep the header of the first file only, then append every data row.
    """
    awk 'FNR > 1 || NR == 1' cells/*/* > endomorph_cells_merged.csv
    awk 'FNR > 1 || NR == 1' summaries/*/* > endomorph_summary_merged.csv
    """
}
