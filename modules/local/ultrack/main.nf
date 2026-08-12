process ULTRACK {
    tag "$meta.id"
    label 'process_high'
    label 'gpu'
    container "docker.io/royerlab/ultrack:0.6.1-cuda12.4"

    input:
    tuple val(meta), val(labels_zarr), val(foreground_zarr), val(contours_zarr), val(raw_zarr)
    path  config_toml

    output:
    tuple val(meta), val(output_segments_zarr), emit: segments
    tuple val(meta), path(output_tracks_csv),   emit: tracks
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    output_tracks_csv   = "${meta.id}_tracks.csv"
    output_segments_zarr = "${meta.id}_segments.zarr"
    def working_dir = "${meta.id}_ultrack_db"
    def config_arg     = config_toml ? "--config-toml ${config_toml}" : ""
    def labels_arg      = labels_zarr    ? "--labels-zarr \"${labels_zarr}\""       : ""
    def foreground_arg  = foreground_zarr ? "--foreground-zarr \"${foreground_zarr}\"" : ""
    def contours_arg    = contours_zarr   ? "--contours-zarr \"${contours_zarr}\""     : ""
    // Only used when contours_zarr is not given -- derive contours from raw
    // intensity data via robust_invert instead
    def derive_contours_arg = (!contours_zarr && raw_zarr) ? "--derive-contours-from-raw \"${raw_zarr}\"" : ""
    def foreground_channel_arg = meta.foreground_channel != null ? "--foreground-channel ${meta.foreground_channel}" : ""
    def contours_channel_arg   = meta.contours_channel   != null ? "--contours-channel ${meta.contours_channel}"     : ""
    def raw_channel_arg        = meta.raw_channel        != null ? "--raw-channel ${meta.raw_channel}"               : ""
    def scale_z_arg = meta.scale_z != null ? "--scale-z ${meta.scale_z}" : ""
    def scale_y_arg = meta.scale_y != null ? "--scale-y ${meta.scale_y}" : ""
    def scale_x_arg = meta.scale_x != null ? "--scale-x ${meta.scale_x}" : ""
    """
    run_ultrack.py \\
        ${labels_arg} \\
        ${foreground_arg} \\
        ${contours_arg} \\
        ${derive_contours_arg} \\
        ${foreground_channel_arg} \\
        ${contours_channel_arg} \\
        ${raw_channel_arg} \\
        ${config_arg} \\
        --working-dir ${working_dir} \\
        --output-tracks-csv ${output_tracks_csv} \\
        --output-segments-zarr ${output_segments_zarr} \\
        ${scale_z_arg} \\
        ${scale_y_arg} \\
        ${scale_x_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ultrack: \$(python -c "import ultrack; print(getattr(ultrack, '__version__', '0.6.1'))")
    END_VERSIONS
    """
}
