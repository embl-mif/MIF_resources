process PYMIF_ENDOMORPH_REGIONS {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(labels_zarr)

    output:
    tuple val(meta), val(output_zarr_path), emit: zarr
    path staged_output_zarr, optional: true, emit: output_zarr_dir // to track tiles that are completed and skip them
    path "versions.yml"                       , emit: versions

    script:
    output_zarr = meta.endomorph_regions_output_name ?: "${meta.id}_endomorph_regions.ome.zarr"
    staged_output_zarr = new File(output_zarr.toString()).isAbsolute() ? null : output_zarr
    // Resolved here, not in the output block: on -resume the cached context's `task`
    // has no workDir, so `task.workDir` in an output declaration is null.
    output_zarr_path = task.workDir.resolve(output_zarr.toString()).toString()

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def levelArg = meta.endomorph_labels_resolution_level != null ? "--resolution-level ${meta.endomorph_labels_resolution_level}" : ""
    def timepointArg = meta.endomorph_timepoint != null ? "--timepoint ${meta.endomorph_timepoint}" : ""
    def ringDepthsArg = meta.endomorph_ring_depths ? "--ring-depths ${shellQuote(meta.endomorph_ring_depths)}" : ""
    def borderMarginArg = meta.endomorph_border_margin != null ? "--border-margin ${meta.endomorph_border_margin}" : ""

    """
    if [ -d ${shellQuote(output_zarr)} ] && [ -e ${shellQuote(output_zarr)}/zarr.json -o -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "Output already exists, skipping endo_morph regions: ${output_zarr}"
    else
        mkdir -p \$(dirname ${shellQuote(output_zarr)})

        create_endomorph_regions.py \
            --labels "${labels_zarr}" \
            --output "${output_zarr}" \
            ${levelArg} \
            ${timepointArg} \
            ${ringDepthsArg} \
            ${borderMarginArg}
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
