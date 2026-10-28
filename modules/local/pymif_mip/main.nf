process PYMIF_MIP {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/embl-mif/pymif:0.5.3'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), path(output_mip_link), emit: zarr
    path "versions.yml"                   , emit: versions

    script:
    output_mip = meta.zarr_output ?: "${meta.id}_mip.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_mip_link = new File(output_mip.toString()).name

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    // Keep this script indented with spaces only: Nextflow dedents it, and stray tabs
    // leave the END_VERSIONS heredoc terminator indented, so bash never closes it.
    """
    if [ -d ${shellQuote(output_mip)} ] && [ -e ${shellQuote(output_mip)}/zarr.json -o -e ${shellQuote(output_mip)}/.zgroup ]; then
        echo "Output already exists, skipping mip: ${output_mip}"
    else
        mkdir -p \$(dirname ${shellQuote(output_mip)})

        create_pymif_mip.py \
            --input "${input_zarr}" \
            --output "${output_mip}"
    fi

    if [ ! -e ${shellQuote(output_mip)}/zarr.json ] && [ ! -e ${shellQuote(output_mip)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_mip}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_mip)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
