process PYMIF_MIP {
    tag "${meta.id}"
    label 'process_medium'
    container 'ghcr.io/grinic/pymif:2026.7.4'

    input:
    tuple val(meta), path(input_zarr)

    output:
    tuple val(meta), val("${task.workDir.resolve(output_mip.toString())}"), emit: zarr
    path staged_output_mip, optional: true, emit: output_dir // to track tiles that are completed and skip them
    path "versions.yml"                       , emit: versions

    script:
	output_mip = meta.mip_output_name ?: "${meta.id}_mip.ome.zarr"
    staged_output_mip = new File(output_mip.toString()).isAbsolute() ? null : output_mip
    
	def shellQuote = { value ->
        "'${value.toString().replace("'", "'\"'\"'")}'"
    }

	"""	
    if [ -d ${shellQuote(output_mip)} ] && [ -e ${shellQuote(output_mip)}/zarr.json -o -e ${shellQuote(output_mip)}/.zgroup ]; then
         echo "Output already exists, skipping mip: ${output_mip}"
    else
        mkdir -p \$(dirname ${shellQuote(output_mip)})
        
		create_pymif_mip.py \
			--input "${input_zarr}" \
			--output "${output_mip}"
    fi
 
 	if [ ! -e ${shellQuote(output_mip)}/zarr.json ] && [ ! -e ${shellQuote(output_mip)}/.zgroup ]; then
         echo "ERROR: expected zarr store not found: ${output_mip}" >&2
         exit 1
    fi  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymif: \$(python -c 'import pymif; print(pymif.__version__)' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
