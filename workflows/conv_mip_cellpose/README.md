# conv_mip_cellpose

Runs each dataset through three stages:

1. Convert to OME-Zarr with `zarrsmith convert`.
2. Build a multiscale 2D maximum-intensity projection (MIP) over Z.
3. Segment the MIP with Cellpose-SAM.

## Input CSV

One dataset per row. Required column:

- `input`: path to the input microscopy dataset.

### Output & Cellpose columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Name of the converted OME-Zarr. An absolute path is kept as-is; anything else is used as a plain filename and stays in the Nextflow work directory. |
| `mip_output` | Name of the MIP OME-Zarr, same rule as `output`. Defaults to `<dataset_id>_mip.ome.zarr`. |
| `labels_output` | Name of the Cellpose label OME-Zarr, same rule as `output`. Defaults to `<dataset_id>_labels.ome.zarr`. |
| `cellpose_diameter` | Expected object diameter in pixels (integer). |
| `cellpose_niter` | Number of Cellpose iterations (integer). |
| `cellpose_channels` | Channels used for segmentation. |
| `cellpose_halo` | Segmentation-context halo in pixels. Defaults to `20`. |
| `resolution_level` | Pyramid level used for segmentation. Defaults to `--resolution_level` (0). |
| `do_3D` | Any non-empty value enables 3D segmentation; leave blank to disable. |

### Conversion columns (`zarrsmith convert`)

Same set as [conv_cellpose](../conv_cellpose/README.md#conversion-columns-zarrsmith-convert):
`microscope`, `chunk_size`, `max_size(MB)`, `scene_index`, `zarr_format`,
`downscale_factor`, `channel_colors`, `channel_names`, `num_levels`, `subset`,
`shards`, `shard_target_mb`, `shard_exclude_axes`. `chunk_size` takes
precedence over `max_size(MB)` when both are set. Paths containing commas must
be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,cellpose_diameter,cellpose_channels
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,30,0
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,25,0
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--resolution_level` | `0` | Pyramid level when a row has none. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/nf-modules/workflows/conv_mip_cellpose/main.nf \
    -lib /path/to/nf-modules/lib \
    --input_csv /path/to/input.csv \
    -resume
```

For batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit:

```bash
mkdir -p logs
sbatch template_run.sh
```
