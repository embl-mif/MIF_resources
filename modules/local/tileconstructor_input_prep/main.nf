process TILECONSTRUCTOR_INPUT_PREP {
    tag "generating input for tileconstructor"
    label 'process_low'
    conda "${moduleDir}/environment.yml"


    input:
    val input_dir
    val resolution_level

    output:
    path "tile_constructor_config.csv",    emit: csv
    path "versions.yml",                   emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    script:
    """
    tileconstructor_config.py \
        --input-dir ${input_dir} \
        --output-csv tile_constructor_config.csv \
        --resolution-level ${resolution_level} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
