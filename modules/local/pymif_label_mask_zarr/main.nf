process PYMIF_LABEL_MASK_ZARR {
    tag "${meta.id}"
    label 'process_high'
    container 'ghcr.io/grinic/pymif:2026.9.1'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), val("${task.workDir.resolve(output_mask.toString())}"), emit: zarr
    path staged_output_mask, optional: true, emit: output_dir // to track tiles that are completed and skip them
    path staged_regionprops_csv, optional: true, emit: regionprops
    path "versions.yml"                       , emit: versions

    script:
    output_mask = meta.mask_output_name ?: "${meta.id}_masked.ome.zarr"
    staged_output_mask = new File(output_mask.toString()).isAbsolute() ? null : output_mask
    // create_pymif_label_mask.py derives the measurements filename from
    // Path(output).stem, which (like pathlib for any "*.zarr" path) strips
    // only the trailing ".zarr" suffix -- not ".ome.zarr" -- so mirror that
    // here rather than assuming a fixed "*.ome.zarr" naming convention.
    def output_mask_stem = output_mask.toString().replaceAll(/\.zarr$/, '')
    regionprops_csv = "${output_mask_stem}_regionprops.csv"
    staged_regionprops_csv = staged_output_mask ? regionprops_csv : null

    def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

    def otsuMultiplierArg = meta.mask_otsu_multiplier != null ? "--otsu-multiplier ${meta.mask_otsu_multiplier}" : ""
    def minAreaArg = meta.mask_min_area != null ? "--min-area ${meta.mask_min_area}" : ""
    def maxAreaArg = meta.mask_max_area != null ? "--max-area ${meta.mask_max_area}" : ""
    def intensityChannelsArg = meta.mask_intensity_channels ? "--intensity-channels ${shellQuote(meta.mask_intensity_channels)}" : ""
    def filtersArg = meta.mask_filters ? "--filters ${shellQuote(meta.mask_filters)}" : ""
    def resolutionLevelArg = meta.mask_resolution_level != null ? "--resolution-level ${meta.mask_resolution_level}" : ""

    """
    if [ -d ${shellQuote(output_mask)} ] && [ -e ${shellQuote(output_mask)}/zarr.json -o -e ${shellQuote(output_mask)}/.zgroup ]; then
        echo "Output already exists, skipping label mask: ${output_mask}"
    else
        mkdir -p \$(dirname ${shellQuote(output_mask)})

        create_pymif_label_mask.py \
            --input "${input_zarr}" \
            --output "${output_mask}" \
            --mask-channel ${meta.mask_channel} \
            ${resolutionLevelArg} \
            ${otsuMultiplierArg} \
            ${minAreaArg} \
            ${maxAreaArg} \
            ${intensityChannelsArg} \
            ${filtersArg}
    fi

    if [ ! -e ${shellQuote(output_mask)}/zarr.json ] && [ ! -e ${shellQuote(output_mask)}/.zgroup ]; then
        echo "ERROR: expected zarr store not found: ${output_mask}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
