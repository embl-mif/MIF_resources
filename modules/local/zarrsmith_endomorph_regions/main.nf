process ZARRSMITH_ENDOMORPH_REGIONS {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/embl-mif/zarrsmith:0.6.1'

    input:
    // nucleus_zarr: the clean-channel zarr, whose "<channel>_clean" nuclei are never gap.
    tuple val(meta), path(labels_zarr), path(intensity_zarr), path(nucleus_zarr, stageAs: 'nucleus/*')

    output:
    tuple val(meta), path(output_zarr_link), emit: zarr
    tuple val(meta), path(empty_qc_json)   , emit: qc
    path "versions.yml"                    , emit: versions

    script:
    output_zarr = meta.zarr_output ?: "${meta.id}_endomorph_regions.ome.zarr"
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
    def emptyNucleusChannelArg = meta.endomorph_empty_nucleus_channel ? "--empty-nucleus-channel ${shellQuote(meta.endomorph_empty_nucleus_channel)}" : ""
    def emptyNucleusFactorArg = meta.endomorph_empty_nucleus_factor != null ? "--empty-nucleus-factor ${meta.endomorph_empty_nucleus_factor}" : ""
    def emptyCoreMinAreaArg = meta.endomorph_empty_core_min_area != null ? "--empty-core-min-area ${meta.endomorph_empty_core_min_area}" : ""
    def emptyCellPercentileArg = meta.endomorph_empty_cell_percentile != null ? "--empty-cell-percentile ${meta.endomorph_empty_cell_percentile}" : ""
    def emptyCoreThresholdArg = meta.endomorph_empty_core_threshold != null ? "--empty-core-threshold ${meta.endomorph_empty_core_threshold}" : ""
    def emptyEdgeThresholdArg = meta.endomorph_empty_edge_threshold != null ? "--empty-edge-threshold ${meta.endomorph_empty_edge_threshold}" : ""
    def emptyMaxFloorRatioArg = meta.endomorph_empty_max_floor_ratio != null ? "--empty-max-floor-ratio ${meta.endomorph_empty_max_floor_ratio}" : ""

    """
    mkdir -p \$(dirname ${shellQuote(output_zarr)})

    create_endomorph_regions.py \
        --labels "${labels_zarr}" \
        --intensity "${intensity_zarr}" \
        --empty-nucleus-zarr "${nucleus_zarr}" \
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
        ${emptyMinContrastArg} \
        ${emptyNucleusChannelArg} \
        ${emptyNucleusFactorArg} \
        ${emptyCoreMinAreaArg} \
        ${emptyCellPercentileArg} \
        ${emptyCoreThresholdArg} \
        ${emptyEdgeThresholdArg} \
        ${emptyMaxFloorRatioArg}

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
        zarrsmith: \$(python -c 'from importlib.metadata import version; print(version("zarrsmith"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
