process PYMIF_OTSU_LABEL {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/embl-mif/pymif:0.5.0'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), path(output_labels_link), emit: zarr
    path "versions.yml"                      , emit: versions

    script:
    output_labels = meta.otsu_output_name ?: "${meta.id}_otsu_labels.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_labels_link = new File(output_labels.toString()).name

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

    ${Utils.shellLinkOutput(output_labels)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
