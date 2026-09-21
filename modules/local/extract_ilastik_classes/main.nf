process EXTRACT_ILASTIK_CLASSES {
    tag "$project"
    label 'process_low'
    container "docker.io/biocontainers/ilastik:1.4.2_cv1"

    input:
    path project

    output:
    path "ilp_classes.json",  emit: classes
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ilastik_python = "/opt/ilastik-1.4.2-Linux/bin/python3"
    """
    ${ilastik_python} ${moduleDir}/resources/usr/bin/extract_ilastik_classes.py \\
        --project "${project}" > ilp_classes.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        h5py: \$(${ilastik_python} -c "import h5py; print(h5py.__version__)")
    END_VERSIONS
    """
}
