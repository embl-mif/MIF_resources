process REBUILD_PYRAMID {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(output_zarr)

    output:
    tuple val(meta), val(output_zarr), emit: zarr
    path "versions.yml",                emit: versions

    script:
    def n_threads = task.cpus ?: 1
    """
    rebuild_pyramid.py \
        --output-zarr "${output_zarr}" \
        --downsample-method Average \
        --n-threads ${n_threads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squirrel: \$(python -c "import squirrel; print(getattr(squirrel, '__version__', 'unknown'))" 2>/dev/null || echo "git-installed")
        zarr: \$(python -c "import zarr; print(zarr.__version__)")
    END_VERSIONS
    """
}
