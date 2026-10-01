process RUN_ULTRACK {
    tag "$meta.id"
    label 'process_high'
    label 'gpu'
    container "docker.io/royerlab/ultrack:0.6.1-cuda12.4"

    input:
    tuple val(meta),
          path(local_labels_zarr),     val(labels_zarr_uri),
          path(local_foreground_zarr), val(foreground_zarr_uri),
          path(local_contours_zarr),   val(contours_zarr_uri),
          path(local_raw_zarr),        val(raw_zarr_uri),
          path(output_tracked_zarr)
    path  config_toml

    output:
    tuple val(meta), path(output_tracked_zarr), emit: tracked_zarr
    tuple val(meta), path(output_tracks_csv),   emit: tracks
    tuple val(meta), path(output_graph),        emit: graph
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    output_tracks_csv = "${meta.id}_tracks.csv"
    output_graph      = "${meta.id}_graph.npy"

    def working_dir = "${meta.id}_ultrack_db"
    def n_threads = task.cpus ?: 1

    // Each of labels/foreground/contours/raw is independently optional -- a
    // caller that doesn't have one passes the NO_FILE sentinel as its local
    // path and null/'' as its uri (same convention as CELLPOSE_SAM_ZARR_MIF's
    // model input). Resolving back to plain labels_zarr/foreground_zarr/etc.
    // here keeps every arg line below completely unchanged.
    def resolveZarr = { localPath, uri -> uri ?: (localPath.name != 'NO_FILE' ? localPath.toString() : null) }
    def labels_zarr     = resolveZarr(local_labels_zarr, labels_zarr_uri)
    def foreground_zarr = resolveZarr(local_foreground_zarr, foreground_zarr_uri)
    def contours_zarr   = resolveZarr(local_contours_zarr, contours_zarr_uri)
    def raw_zarr         = resolveZarr(local_raw_zarr, raw_zarr_uri)

    def t_min_arg = meta.t_min != null ? "--t-min ${meta.t_min}" : ""
    def t_max_arg = meta.t_max != null ? "--t-max ${meta.t_max}" : ""
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
    def min_area_arg      = meta.min_area      != null ? "--min-area ${meta.min_area}"           : ""
    def max_area_arg      = meta.max_area      != null ? "--max-area ${meta.max_area}"            : ""
    def max_distance_arg  = meta.max_distance  != null ? "--max-distance ${meta.max_distance}"    : ""
    def distance_weight_arg = meta.distance_weight != null ? "--distance-weight ${meta.distance_weight}" : ""
    def solution_gap_arg  = meta.solution_gap  != null ? "--solution-gap ${meta.solution_gap}"     : ""
    def time_limit_arg    = meta.time_limit    != null ? "--time-limit ${meta.time_limit}"         : ""
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
        ${t_min_arg} \\
        ${t_max_arg} \\
        ${config_arg} \\
        --working-dir ${working_dir} \\
        --output-tracks-csv ${output_tracks_csv} \\
        --output-graph ${output_graph} \\
        --output-segments-zarr ${output_tracked_zarr} \\
        --n-threads ${n_threads} \\
        ${min_area_arg} \\
        ${max_area_arg} \\
        ${max_distance_arg} \\
        ${distance_weight_arg} \\
        ${solution_gap_arg} \\
        ${time_limit_arg} \\
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
