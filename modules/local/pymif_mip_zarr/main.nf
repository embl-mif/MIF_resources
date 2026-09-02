process PYMIF_MIP_ZARR {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.7.4'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), val(meta.mip_output_name), emit: zarr
    path "versions.yml"                       , emit: versions

    script:
    """
    create_pymif_mip.py \
        --input "${input_zarr}" \
        --output "${meta.mip_output_name}"
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}