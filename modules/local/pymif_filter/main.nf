process PYMIF_FILTER {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(labels_zarr), path(csv_file)

    output:
    tuple val(meta), val("${task.workDir.resolve(output_labels.toString())}"), val("${task.workDir.resolve(output_csv.toString())}"), emit: zarr_csv
    path staged_output_labels, optional: true, emit: output_zarr // to track tiles that are completed and skip them
    path staged_output_csv, optional: true, emit: output_csv // to track tiles that are completed and skip them
    path "versions.yml"                       , emit: versions

    script:
    output_labels = meta.filter_output_name ?: "${meta.id}_filtered_labels.ome.zarr"
    staged_output_labels = new File(output_labels.toString()).isAbsolute() ? null : output_labels
    output_csv = meta.filter_csv_output_name ?: "${meta.id}_filtered_labels.csv"
    staged_output_csv = new File(output_csv.toString()).isAbsolute() ? null : output_csv

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def minAreaArg = meta.filter_min_area != null ? "--min-area ${meta.filter_min_area}" : ""
    def maxAreaArg = meta.filter_max_area != null ? "--max-area ${meta.filter_max_area}" : ""
    def filtersArg = meta.filter_properties ? "--filters ${shellQuote(meta.filter_properties)}" : ""

    """
    if [ -d ${shellQuote(output_labels)} ] && [ -e ${shellQuote(output_labels)}/zarr.json -o -e ${shellQuote(output_labels)}/.zgroup ]; then
        echo "Output already exists, skipping label filtering: ${output_labels}"
    else
        mkdir -p \$(dirname ${shellQuote(output_labels)})
        mkdir -p \$(dirname ${shellQuote(output_csv)})

        filter_pymif_labels.py \
            --labels "${labels_zarr}" \
            --csv "${csv_file}" \
            --output "${output_labels}" \
            --output-csv "${output_csv}" \
            ${minAreaArg} \
            ${maxAreaArg} \
            ${filtersArg}
    fi

    if [ ! -e ${shellQuote(output_labels)}/zarr.json ] && [ ! -e ${shellQuote(output_labels)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_labels}" >&2
        exit 1
    fi

    if [ ! -e ${shellQuote(output_csv)} ]; then
        echo "ERROR: expected filtered measurements CSV not found: ${output_csv}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
