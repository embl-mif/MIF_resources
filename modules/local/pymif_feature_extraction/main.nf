process PYMIF_FEATURE_EXTRACTION {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(intensity_zarr), path(labels_zarr)

    output:
    tuple val(meta), val("${task.workDir.resolve(output_csv.toString())}"), emit: csv
    path staged_output_csv, optional: true, emit: output_file // to track tiles that are completed and skip them
    path "versions.yml"                    , emit: versions

    script:
    output_csv = meta.features_output_name ?: "${meta.id}_features.csv"
    staged_output_csv = new File(output_csv.toString()).isAbsolute() ? null : output_csv

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
