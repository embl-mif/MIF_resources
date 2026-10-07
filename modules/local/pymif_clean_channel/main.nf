process PYMIF_CLEAN_CHANNEL {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.10.0'

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    path "versions.yml"                    , emit: versions

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}_clean.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_zarr_link = new File(output_zarr.toString()).name

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    """
    if [ -d ${shellQuote(output_zarr)} ] && [ -e ${shellQuote(output_zarr)}/zarr.json -o -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "Output already exists, skipping clean channel: ${output_zarr}"
    else
        mkdir -p \$(dirname ${shellQuote(output_zarr)})

        create_pymif_clean_channel.py \
            --intensity "${intensity_zarr}" \
            --labels "${labels_zarr}" \
            --output "${output_zarr}" \
            --channel ${meta.clean_channel}
    fi

    if [ ! -e ${shellQuote(output_zarr)}/zarr.json ] && [ ! -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_zarr}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_zarr)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
