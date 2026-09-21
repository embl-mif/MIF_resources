process SAMPLE_INPUT_PREP {
    tag "Prepare sample input for ${exp_dir}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path exp_dir
    val zarr_ver
    val sample_table

    output:
    path "${sample_table}", emit: csv
    path "versions.yml",     emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    generate_sample_table.py \
        --exp ${exp_dir} \
        --zarr-version ${zarr_ver} \
        --out ${sample_table} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
