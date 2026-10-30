process ZARRSMITH_CLEAN_CHANNEL {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/embl-mif/zarrsmith:0.6.1'

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
    out=${shellQuote(output_zarr)}
    if [ -d "\$out" ] && { grep -q '"ome"' "\$out/zarr.json" 2>/dev/null || grep -q '"multiscales"' "\$out/.zattrs" 2>/dev/null; }; then
        echo "Output already exists, skipping clean channel: ${output_zarr}"
    else
        if [ -e "\$out" ]; then
            echo "Incomplete output found (no OME metadata), removing: ${output_zarr}"
            rm -rf "\$out"
        fi
        mkdir -p "\$(dirname "\$out")"

        create_zarrsmith_clean_channel.py \
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
        zarrsmith: \$(python -c 'from importlib.metadata import version; print(version("zarrsmith"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
