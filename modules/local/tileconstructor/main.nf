process TILE_CONSTRUCTOR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml" 

    input:
    tuple val(meta), val(image)

    output:
    tuple val(meta), path("*.csv"), emit: tiles

    script:
    def output_csv = "${meta.id}_tiles.csv"
    def resolution_level_arg = meta.resolution_level != null ? "--resolution-level ${meta.resolution_level}" : ""
    def tile_size_x_arg = meta.tile_size_x != null ? "--tile-size-x ${meta.tile_size_x}" : ""
    def tile_size_y_arg = meta.tile_size_y != null ? "--tile-size-y ${meta.tile_size_y}" : ""
    def tile_size_z_arg = meta.tile_size_z != null ? "--tile-size-z ${meta.tile_size_z}" : ""
    def tile_overlap_arg = meta.tile_overlap != null ? "--tile-overlap ${meta.tile_overlap}" : ""
    def x_min_arg = meta.x_min != null ? "--x-min ${meta.x_min}" : ""
    def x_max_arg = meta.x_max != null ? "--x-max ${meta.x_max}" : ""
    def y_min_arg = meta.y_min != null ? "--y-min ${meta.y_min}" : ""
    def y_max_arg = meta.y_max != null ? "--y-max ${meta.y_max}" : ""
    def z_min_arg = meta.z_min != null ? "--z-min ${meta.z_min}" : ""
    def z_max_arg = meta.z_max != null ? "--z-max ${meta.z_max}" : ""
    def use_physical_units_arg = meta.use_physical_units ? "--use-physical-units" : ""

    if (meta.tile_size_x == null || meta.tile_size_y == null) {
        throw new IllegalArgumentException("tile_size_x and tile_size_y are required in meta for TILE_CONSTRUCTOR")
    }

    """
    tile_constructor_2.py \\
        --source-dataset-id "${meta.id}" \\
        --input-image-path \"${image}\" \\
        --output-csv ${output_csv} \\
        ${tile_size_x_arg} \\
        ${tile_size_y_arg} \\
        ${tile_size_z_arg} \\
        ${tile_overlap_arg} \\
        ${resolution_level_arg} \\
        ${x_min_arg} \\
        ${x_max_arg} \\
        ${y_min_arg} \\
        ${y_max_arg} \\
        ${z_min_arg} \\
        ${z_max_arg} \\
        ${use_physical_units_arg}
    """

    stub:
    def output_csv = "${meta.id}_tiles.csv"
    """
    cat > ${output_csv} <<-END_TILE_CSV
    dataset_id,source_dataset_id,tile_id,input_uri,x_min,x_max,y_min,y_max,z_min,z_max,resolution_level
    ${meta.id}__tile_000001,${meta.id},tile_000001,${image},0,1,0,1,,,0
    END_TILE_CSV
    """
}

