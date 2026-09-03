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
    container "docker://registry.git.embl.org/grp-cba/containers/cellposesam-zarr:4.1.1"

    input:
    tuple val(meta), val(image_uri)

    output:
    tuple val(meta), val(meta.labels_output_name), emit: labels

    script:
    def image = image_uri ?: local_image
    // pretrained_cellpose_model is always a staged, non-null Path (even the NO_MODEL_FILE
    // sentinel from assets/NO_MODEL_FILE when no custom model was requested -- see
    // tiled_cellpose_sam_zarr subworkflow), so presence must be checked by name,
    // not truthiness (same pitfall documented in locallabelmatcher/main.nf). It's a
    // distinct filename from the local_image sentinel (assets/NO_FILE) on purpose --
    // two `path` inputs resolving to the same filename in one task otherwise trigger a
    // Nextflow "input file name collision" error.
    def cellpose_model_arg = meta.pretrained_cellpose_model != null ? "-m ${meta.pretrained_cellpose_model}" : ""

    def x_min_arg = Utils.optionalCliArg("--x-min", meta.x_min)
    def x_max_arg = Utils.optionalCliArg("--x-max", meta.x_max)
    def y_min_arg = Utils.optionalCliArg("--y-min", meta.y_min)
    def y_max_arg = Utils.optionalCliArg("--y-max", meta.y_max)
    def z_min_arg = Utils.optionalCliArg("--z-min", meta.z_min)
    def z_max_arg = Utils.optionalCliArg("--z-max", meta.z_max)
    def channels_arg = Utils.optionalCliArg("--segmentation_channels", meta.segmentation_channels)
    def halo_arg = Utils.optionalCliArg("--halo", meta.halo)
    def resolution_level_arg = Utils.optionalCliArg("--resolution-level", meta.resolution_level)
    def use_physical_units_arg = meta.use_physical_units ? "--use-physical-units" : ""
    
    def diameter_arg = meta.diameter != null ? "--diameter ${meta.diameter}" : ""
    def niter_arg = meta.niter != null ? "--niter ${meta.niter}" : ""
    def flow_threshold_arg = meta.flow_threshold != null ? "--flow-threshold ${meta.flow_threshold}" : ""
    def cellprob_threshold_arg = meta.cellprob_threshold != null ? "--cellprob-threshold ${meta.cellprob_threshold}" : ""

    def args = task.ext.args ?: ''

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
        ${cellpose_model_arg} \
        ${diameter_arg} \
        ${niter_arg} \
        ${flow_threshold_arg} \
        ${cellprob_threshold_arg} \
        ${args}
    """

    stub:
    def output_zarr = "${meta.id}_labels.ome.zarr"
    """
    mkdir -p ${output_zarr}
    touch ${output_zarr}/.zgroup
    """

}
