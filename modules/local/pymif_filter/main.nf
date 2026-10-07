process PYMIF_FILTER {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.10.0'

    input:
    tuple val(meta), path(labels_zarr), path(csv_file)

    output:
    tuple val(meta), path(output_labels_link), path(output_csv_link), emit: zarr_csv
    path "versions.yml"                                             , emit: versions

    script:
    output_labels = meta.zarr_output ?: "${meta.id}_filtered_labels.ome.zarr"
    output_csv = meta.csv_output ?: "${meta.id}_filtered_labels.csv"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_labels_link = new File(output_labels.toString()).name
    output_csv_link = new File(output_csv.toString()).name

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

    ${Utils.shellLinkOutput(output_labels)}
    ${Utils.shellLinkOutput(output_csv)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
