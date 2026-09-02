# conv_mip_cellpose_edges

`conv_mip_cellpose_edges` runs one microscopy dataset per input CSV row through
three stages:

1. Convert the source dataset to OME-Zarr with `pymif 2zarr`.
2. Create a multiscale, two-dimensional maximum-intensity-projection (MIP)
   OME-Zarr by projecting the Z axis.
3. Run Cellpose-SAM on the MIP and write an OME-Zarr label image.

## Input CSV

The CSV must contain one dataset per row. The required columns are:

- `input`: path to the input microscopy dataset.
- `output`: name of the converted OME-Zarr directory.

The workflow also reads these optional columns:

| CSV column | Purpose |
| --- | --- |
| `mip_output` | Name of the MIP ome.zarr dataset. Default to '${dataset_id}_mip.ome.zarr' |
| `dataset_id` | Dataset identifier used in task names and output names. Defaults to the input dataset basename. |
| `cellpose_diameter` | Cellpose object diameter in pixels. The current workflow parses it as an integer. |
| `cellpose_niter` | Number of Cellpose iterations. Must be an integer. |
| `cellpose_channels` | Channels to be used for cellpose segmentation. |

The remaining optional conversion columns are passed to `pymif 2zarr`:

| CSV column | PyMIF argument |
| --- | --- |
| `microscope` | `--microscope` |
| `chunk_size` | `--chunk_size` |
| `max_size(MB)` | `--max_size` |
| `scene_index` | `--scene_index` |
| `zarr_format` | `--zarr_format` |
| `downscale_factor` | `--downscale_factor` |
| `channel_colors` | `--channel_colors` |
| `channel_names` | `--channel_names` |
| `num_levels` | `--num_levels` |
| `subset` | `--subset` |

When both `chunk_size` and `max_size(MB)` are present, `chunk_size` takes
precedence.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,max_size(MB),scene_index,zarr_format,downscale_factor,channel_colors,channel_names,num_levels,subset,cellpose_diameter,cellpose_niter
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,,0,2,2,,"DAPI,GFP",4,,30,200
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,100,0,2,2,,"DAPI,GFP",4,,25,150
```

Empty optional fields are passed through as unset values. Paths containing
commas must be CSV-quoted.

## Output locations

Converted datasets are written to:

```text
<outdir>/zarr/<output basename>
```

The MIP and Cellpose modules emit their OME-Zarr directories as Nextflow
channels. Unless a `publishDir` rule is configured for those processes, their
artifacts remain in the Nextflow work directory.

## Running the workflow

The workflow configuration is set up for the EMBL Heidelberg Slurm cluster and
uses Apptainer for its containers.

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_mip_cellpose_edges/main.nf \
    --input_csv /path/to/input.csv \
    --outdir /scratch/$USER/conv_mip_cellpose_edges/results \
    -resume
```

For batch submission, copy `template_run.sh`, update the marked paths and email
address, create its `logs` directory, and submit it:

```bash
mkdir -p logs
sbatch template_run.sh
```

## Module contracts

The workflow expects the local modules to exchange these tuples:

```groovy
// PYMIF_CONVERSION.out.zarr
tuple(meta, converted_zarr)

// MIP_ZARR.out.zarr
tuple(meta, mip_zarr)

// CELLPOSE_SAM_ZARR_WIP input
tuple(meta, mip_zarr)
```

`modules/local/mip_zarr/main.nf` must therefore exist and expose an output
named `zarr` before this workflow can run.
