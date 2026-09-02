process CELLPOSE_SAM_ZARR_MIF {
    tag "$meta.id"
    label "process_high"
    label "gpu"

    // Catch GPU OOM errors
    // FIXME: Also add retries for out of time
    errorStrategy {
        def errFile = new File("${task.workDir}/.command.err")
        def stderr = errFile.exists() ? errFile.text : ""
        def lowered = stderr.toLowerCase()
        def isGPUOOM = [
            'out of memory',
            'cuda_error_out_of_memory',
            'cudaerrormemoryallocation',
            'resourceexhaustederror',
            'cuda out of memory',
            'cudnn fails to initialize',
            'failed to allocate',
        ].any { pattern -> lowered.contains(pattern) }
        def isGPUInitFailure = [
            'gpu_init_error',
            'cuda unknown error',
            'cuda initialization',
            'no cuda gpus are available',
            'torch not compiled with cuda',
            'cuda driver version is insufficient',
            'neither torch cuda nor mps version not installed/working',
        ].any { pattern -> lowered.contains(pattern) }

        if (isGPUOOM || isGPUInitFailure) {
            def retryCause = isGPUOOM ? 'GPU_OOM' : 'GPU_INIT_FAILURE'
            def lines = stderr ? stderr.readLines() : []
            def interesting = lines.findAll { line ->
                def l = line.toLowerCase()
                l.contains('gpu_init_error') ||
                l.contains('cuda') ||
                l.contains('torch') ||
                l.contains('cellpose') ||
                l.contains('out of memory') ||
                l.contains('traceback')
            }
            def excerpt = (interesting ? interesting : lines).takeRight(25).join('\n')
            log.warn(
                "Retrying ${task.process} [${task.tag}] attempt ${task.attempt}/${task.maxRetries + 1} " +
                "(cause=${retryCause}, jobId=${task.native_id}, workDir=${task.workDir})\n" +
                "--- .command.err excerpt ---\n${excerpt}\n--- end excerpt ---"
            )
            return 'retry'
        }

        return 'terminate'
    }
    maxRetries 2
    
    // Set up environment
    // FIXME: Replace conda with a container
    container "docker://registry.git.embl.org/grp-cba/containers/cellposesam-zarr:4.1.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), val(meta.labels_output_name), emit: labels

    script:
    def x_min_arg = meta.x_min != null ? "--x-min ${meta.x_min}" : ""
    def x_max_arg = meta.x_max != null ? "--x-max ${meta.x_max}" : ""
    def y_min_arg = meta.y_min != null ? "--y-min ${meta.y_min}" : ""
    def y_max_arg = meta.y_max != null ? "--y-max ${meta.y_max}" : ""
    def z_min_arg = meta.z_min != null ? "--z-min ${meta.z_min}" : ""
    def z_max_arg = meta.z_max != null ? "--z-max ${meta.z_max}" : ""
    def channels_arg = meta.segmentation_channels != null ? "--segmentation_channels ${meta.segmentation_channels}" : ""
    def halo_arg = meta.halo != null ? "--halo ${meta.halo}" : ""
    def resolution_level_arg = meta.resolution_level != null ? "--resolution-level ${meta.resolution_level}" : ""
    def use_physical_units_arg = meta.use_physical_units ? "--use-physical-units" : ""
    def diameter = meta.diameter != null ? "--diameter ${meta.diameter}" : ""
    def niter = meta.niter != null ? "--niter ${meta.niter}" : ""
    """
    cellpose_sam_zarr.py \
        --input_zarr "${image}" \
        --output_zarr ${meta.labels_output_name} \
        ${channels_arg} \
        ${x_min_arg} \
        ${x_max_arg} \
        ${y_min_arg} \
        ${y_max_arg} \
        ${z_min_arg} \
        ${z_max_arg} \
        ${halo_arg} \
        ${resolution_level_arg} \
        ${use_physical_units_arg} \
        ${diameter} \
        ${niter}
    """

    stub:
    """
    mkdir -p ${meta.labels_output_name}
    touch ${meta.labels_output_name}/.zgroup
    """

}

