process PYMIF_ENDOMORPH_FEATURES {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr), path(regions_zarr)

    output:
    tuple val(meta), path(output_csv_link), path(output_summary_csv_link), emit: csv
    path "versions.yml"                                                  , emit: versions

    script:
    output_csv = meta.endomorph_cells_output_name ?: "${meta.id}_endomorph_cells.csv"
    output_summary_csv = meta.endomorph_summary_output_name ?: "${meta.id}_endomorph_summary.csv"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_csv_link = new File(output_csv.toString()).name
    output_summary_csv_link = new File(output_summary_csv.toString()).name

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def intensityLevelArg = meta.endomorph_intensity_resolution_level != null ? "--intensity-resolution-level ${meta.endomorph_intensity_resolution_level}" : ""
    def labelsLevelArg = meta.endomorph_labels_resolution_level != null ? "--labels-resolution-level ${meta.endomorph_labels_resolution_level}" : ""
    def timepointArg = meta.endomorph_timepoint != null ? "--timepoint ${meta.endomorph_timepoint}" : ""
    def vcadArg = meta.endomorph_vcad_channel ? "--vcad-channel ${shellQuote(meta.endomorph_vcad_channel)}" : ""
    def myosinArg = meta.endomorph_myosin_channel ? "--myosin-channel ${shellQuote(meta.endomorph_myosin_channel)}" : ""
    def actinArg = meta.endomorph_actin_channel ? "--actin-channel ${shellQuote(meta.endomorph_actin_channel)}" : ""
    def nucleiArg = meta.endomorph_nuclei_channel ? "--nuclei-channel ${shellQuote(meta.endomorph_nuclei_channel)}" : ""
    def emptyChannelsArg = meta.endomorph_empty_channels ? "--empty-channels ${shellQuote(meta.endomorph_empty_channels)}" : ""
    def emptySmoothingArg = meta.endomorph_empty_smoothing != null ? "--empty-smoothing ${meta.endomorph_empty_smoothing}" : ""
    def emptyMinAreaArg = meta.endomorph_empty_min_area != null ? "--empty-min-area ${meta.endomorph_empty_min_area}" : ""
    def emptySeedPercentileArg = meta.endomorph_empty_seed_percentile != null ? "--empty-seed-percentile ${meta.endomorph_empty_seed_percentile}" : ""
    def emptyMinContrastArg = meta.endomorph_empty_min_contrast != null ? "--empty-min-contrast ${meta.endomorph_empty_min_contrast}" : ""

    """
    if [ -e ${shellQuote(output_csv)} ] && [ -e ${shellQuote(output_summary_csv)} ]; then
        echo "Output already exists, skipping endo_morph feature extraction: ${output_csv}"
    else
        mkdir -p \$(dirname ${shellQuote(output_csv)})
        mkdir -p \$(dirname ${shellQuote(output_summary_csv)})

        extract_endomorph_features.py \
            --intensity "${intensity_zarr}" \
            --labels "${labels_zarr}" \
            --regions "${regions_zarr}" \
            --output "${output_csv}" \
            --output-summary "${output_summary_csv}" \
            --dataset-id ${shellQuote(meta.id)} \
            ${intensityLevelArg} \
            ${labelsLevelArg} \
            ${timepointArg} \
            ${vcadArg} \
            ${myosinArg} \
            ${actinArg} \
            ${nucleiArg} \
            ${emptyChannelsArg} \
            ${emptySmoothingArg} \
            ${emptyMinAreaArg} \
            ${emptySeedPercentileArg} \
            ${emptyMinContrastArg}
    fi

    if [ ! -e ${shellQuote(output_csv)} ]; then
        echo "ERROR: expected cells CSV not found: ${output_csv}" >&2
        exit 1
    fi

    if [ ! -e ${shellQuote(output_summary_csv)} ]; then
        echo "ERROR: expected summary CSV not found: ${output_summary_csv}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_csv)}
    ${Utils.shellLinkOutput(output_summary_csv)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
