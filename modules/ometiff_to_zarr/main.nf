process CONVERT_TO_ZARR {
    tag "${meta.exp_name}-${meta.well}-${meta.field}"
    label 'process_medium'
    container "ghcr.io/grinic/pymif:2026.04.09"

    input:
    tuple path(tiff), val(meta)

    output:
    // stdout
    tuple path("${meta.well}_${meta.field}.zarr"), val(meta), emit: zarr

    script:
    def nl_arg = meta.num_layers?.trim() ?: ""
    def cs_arg = meta.chunk_size?.trim() ?: ""

    """
    echo "${meta.exp_name} - ${meta.well} - ${meta.field}"
    ometiff2zarr.py \
                -i '$tiff' \
                -o '${meta.well}_${meta.field}.zarr' \
                -nl '${nl_arg}' \
                -cs '${cs_arg}'
    """
}
