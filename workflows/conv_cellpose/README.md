# conv_cellpose

`conv_cellpose` converts microscopy datasets to OME-Zarr and runs
Cellpose-SAM segmentation on each dataset.

## Input CSV

The CSV must contain one dataset per row. The required columns are:

- `input`: path to the input microscopy dataset.
- `output`: name of the converted OME-Zarr directory.

The workflow also reads these optional Cellpose columns:

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default label names. Defaults to the input dataset basename. |
| `cellpose_diameter` | Cellpose object diameter in pixels. Parsed as an integer. |
| `cellpose_niter` | Number of Cellpose iterations. Parsed as an integer. |
| `cellpose_channels` | Channels used for Cellpose segmentation. |
| `do_3D` | Set to `true` to enable 3D segmentation. Defaults to `false`. |
| `labels_output` | Name of the output labels OME-Zarr directory. |

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
precedence. Paths containing commas must be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,max_size(MB),scene_index,zarr_format,downscale_factor,channel_colors,channel_names,num_levels,subset,cellpose_diameter,cellpose_niter,cellpose_channels,do_3D,labels_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,,0,2,2,,"DAPI,GFP",4,,30,200,0,False,sample_01_labels.ome.zarr
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,100,0,2,,2,,"DAPI,GFP",4,,25,150,0,False,
```

## Outputs

Converted datasets and labels are written under the configured output directory:

```text
<outdir>/zarr/<output basename>
<outdir>/zarr/<labels_output basename>
```

If `labels_output` is empty, the label directory defaults to
`<dataset_id>_labels.ome.zarr`. Cellpose outputs are emitted through the
Nextflow channel; they remain in the Nextflow work directory unless a
`publishDir` rule is configured.

## Running the workflow

The supplied configuration targets the EMBL Heidelberg Slurm cluster and uses
Apptainer for its containers.

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_cellpose/main.nf \
    --input_csv /path/to/input.csv \
    --outdir /scratch/$USER/conv_cellpose/results \
    -resume
```

For a batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit the script:

```bash
mkdir -p logs
sbatch template_run.sh
```