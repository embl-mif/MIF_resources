process PYMIF_OTSU_LABEL {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), val(output_labels_path), emit: zarr
    path staged_output_labels, optional: true, emit: output_zarr // to track tiles that are completed and skip them
    path "versions.yml"                       , emit: versions

    script:
    output_labels = meta.otsu_output_name ?: "${meta.id}_otsu_labels.ome.zarr"
    staged_output_labels = new File(output_labels.toString()).isAbsolute() ? null : output_labels
    // Resolved here, not in the output block: on -resume the cached context's `task`
    // has no workDir, so `task.workDir` in an output declaration is null.
    output_labels_path = task.workDir.resolve(output_labels.toString()).toString()

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def resolutionLevelArg = meta.otsu_resolution_level != null ? "--resolution-level ${meta.otsu_resolution_level}" : ""
    def otsuMultiplierArg = meta.otsu_multiplier != null ? "--otsu-multiplier ${meta.otsu_multiplier}" : ""

    """
    if [ -d ${shellQuote(output_labels)} ] && [ -e ${shellQuote(output_labels)}/zarr.json -o -e ${shellQuote(output_labels)}/.zgroup ]; then
        echo "Output already exists, skipping Otsu label: ${output_labels}"
    else
        mkdir -p \$(dirname ${shellQuote(output_labels)})

        create_pymif_otsu_label.py \
            --input "${input_zarr}" \
            --output "${output_labels}" \
            --channel ${meta.otsu_channel} \
            ${resolutionLevelArg} \
            ${otsuMultiplierArg}
    fi

    if [ ! -e ${shellQuote(output_labels)}/zarr.json ] && [ ! -e ${shellQuote(output_labels)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_labels}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
