process PREALLOCATE_OUTPUT_ZARR {
    tag "$meta.id"
    label 'process_low'

    // find container
    container 

    input:
    tuple val(meta), val(image)
    val pubdir

    output:
    // Emits the absolute path of the newly created OME-Zarr
    tuple val(meta), val("${pubdir}/${meta.id}_probabilities.ome.zarr"), emit: empty_zarr
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: '--num-channels 2'

    """
    # Create the publication directory on the shared system
    mkdir -p "${pubdir}"

    # Preallocate the skeleton locally in the temporary work directory
    preallocate_zarr.py \\
        --input "${image}" \\
        --output "${meta.id}_probabilities.ome.zarr" \\
        ${args}

    # Publish/copy the directory via bash
    cp -r "${meta.id}_probabilities.ome.zarr" "${pubdir}/"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
    END_VERSIONS
    """
}
