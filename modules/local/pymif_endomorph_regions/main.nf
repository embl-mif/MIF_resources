process PYMIF_ENDOMORPH_REGIONS {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(labels_zarr), path(intensity_zarr)

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    tuple val(meta), path(empty_qc_json)   , emit: qc
    path "versions.yml"                    , emit: versions

    script:
    output_zarr = meta.endomorph_regions_output_name ?: "${meta.id}_endomorph_regions.ome.zarr"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_zarr_link = new File(output_zarr.toString()).name
    // QC diagnostics of the empty-area estimation. Always in the work dir.
    empty_qc_json = "${meta.id}_endomorph_empty_qc.json"

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def levelArg = meta.endomorph_labels_resolution_level != null ? "--resolution-level ${meta.endomorph_labels_resolution_level}" : ""
    def intensityLevelArg = meta.endomorph_intensity_resolution_level != null ? "--intensity-resolution-level ${meta.endomorph_intensity_resolution_level}" : ""
    def timepointArg = meta.endomorph_timepoint != null ? "--timepoint ${meta.endomorph_timepoint}" : ""
    def ringDepthsArg = meta.endomorph_ring_depths ? "--ring-depths ${shellQuote(meta.endomorph_ring_depths)}" : ""
    def borderMarginArg = meta.endomorph_border_margin != null ? "--border-margin ${meta.endomorph_border_margin}" : ""
    def emptyChannelsArg = meta.endomorph_empty_channels ? "--empty-channels ${shellQuote(meta.endomorph_empty_channels)}" : ""
    def emptySmoothingArg = meta.endomorph_empty_smoothing != null ? "--empty-smoothing ${meta.endomorph_empty_smoothing}" : ""
    def emptyMinAreaArg = meta.endomorph_empty_min_area != null ? "--empty-min-area ${meta.endomorph_empty_min_area}" : ""
    def emptySeedPercentileArg = meta.endomorph_empty_seed_percentile != null ? "--empty-seed-percentile ${meta.endomorph_empty_seed_percentile}" : ""
    def emptyMinContrastArg = meta.endomorph_empty_min_contrast != null ? "--empty-min-contrast ${meta.endomorph_empty_min_contrast}" : ""

    """
    mkdir -p \$(dirname ${shellQuote(output_zarr)})

    create_endomorph_regions.py \
        --labels "${labels_zarr}" \
        --intensity "${intensity_zarr}" \
        --output "${output_zarr}" \
        --output-empty-qc "${empty_qc_json}" \
        ${emptyChannelsArg} \
        ${levelArg} \
        ${intensityLevelArg} \
        ${timepointArg} \
        ${ringDepthsArg} \
        ${borderMarginArg} \
        ${emptySmoothingArg} \
        ${emptyMinAreaArg} \
        ${emptySeedPercentileArg} \
        ${emptyMinContrastArg}

    if [ ! -e ${shellQuote(output_zarr)}/zarr.json ] && [ ! -e ${shellQuote(output_zarr)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_zarr}" >&2
        exit 1
    fi

    if [ ! -e ${shellQuote(empty_qc_json)} ]; then
        echo "ERROR: expected empty-area QC file not found: ${empty_qc_json}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_zarr)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
