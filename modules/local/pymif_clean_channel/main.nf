process PYMIF_CLEAN_CHANNEL {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr)

    output:
    tuple val(meta), val("${task.workDir.resolve(output_zarr.toString())}"), emit: zarr
    path staged_output_zarr, optional: true, emit: output_zarr_dir // to track tiles that are completed and skip them
    path "versions.yml"                       , emit: versions

    script:
    output_zarr = meta.clean_output_name ?: "${meta.id}_clean.ome.zarr"
    staged_output_zarr = new File(output_zarr.toString()).isAbsolute() ? null : output_zarr

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
