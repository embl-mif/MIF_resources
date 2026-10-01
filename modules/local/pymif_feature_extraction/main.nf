process PYMIF_FEATURE_EXTRACTION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.10.0'

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr)

    output:
    tuple val(meta), path(output_csv_link), emit: csv
    path "versions.yml"                   , emit: versions

    script:
    output_csv = meta.features_output_name ?: "${meta.id}_features.csv"
    // A `path` output must live in the work dir: an absolute output is exposed there as a
    // symlink named after it (see Utils.shellLinkOutput), a relative one is the output itself.
    output_csv_link = new File(output_csv.toString()).name

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def intensityLevelArg = meta.features_intensity_resolution_level != null ? "--intensity-resolution-level ${meta.features_intensity_resolution_level}" : ""
    def labelsLevelArg = meta.features_labels_resolution_level != null ? "--labels-resolution-level ${meta.features_labels_resolution_level}" : ""
    def channelsArg = meta.features_channels ? "--channels ${shellQuote(meta.features_channels)}" : ""
    def extraPropertiesArg = meta.features_extra_properties ? "--extra-properties ${shellQuote(meta.features_extra_properties)}" : ""
    def extraIntensityPropertiesArg = meta.features_extra_intensity_properties ? "--extra-intensity-properties ${shellQuote(meta.features_extra_intensity_properties)}" : ""

    """
    if [ -e ${shellQuote(output_csv)} ]; then
        echo "Output already exists, skipping feature extraction: ${output_csv}"
    else
        mkdir -p \$(dirname ${shellQuote(output_csv)})

        extract_pymif_features.py \
            --intensity "${intensity_zarr}" \
            --labels "${labels_zarr}" \
            --output "${output_csv}" \
            --dataset-id ${shellQuote(meta.id)} \
            ${intensityLevelArg} \
            ${labelsLevelArg} \
            ${channelsArg} \
            ${extraPropertiesArg} \
            ${extraIntensityPropertiesArg}
    fi

    if [ ! -e ${shellQuote(output_csv)} ]; then
        echo "ERROR: expected features CSV not found: ${output_csv}" >&2
        exit 1
    fi

    ${Utils.shellLinkOutput(output_csv)}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'from importlib.metadata import version; print(version("pymif"))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
