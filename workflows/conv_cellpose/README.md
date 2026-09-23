# conv_cellpose

Converts microscopy datasets to OME-Zarr and segments each with Cellpose-SAM.

## Input CSV

One dataset per row. Required column:

- `input`: path to the input microscopy dataset.

### Output columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Path of the converted OME-Zarr. Must be absolute if set. Defaults to `<default_outdir>/zarr/<dataset_id>.ome.zarr`. |
| `labels_output` | Path of the Cellpose label OME-Zarr. Must be absolute if set. Defaults to `<default_outdir>/labels/<dataset_id>_labels.ome.zarr`. |

### Cellpose columns

| CSV column | Purpose |
| --- | --- |
| `cellpose_diameter` | Expected object diameter in pixels (integer). |
| `cellpose_niter` | Number of Cellpose iterations (integer). |
| `cellpose_channels` | Channels used for segmentation. |
| `cellpose_flow_threshold` | Cellpose flow threshold (float). |
| `cellpose_cellprob_threshold` | Cellpose cell-probability threshold (float). |
| `cellpose_halo` | Segmentation-context halo in pixels. Defaults to `20`. |
| `resolution_level` | Pyramid level used for segmentation. Defaults to `--resolution_level` (0). |
| `use_physical_units` | `true` to interpret sizes/ROI in physical units. Defaults to `false`. |
| `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` | Restrict segmentation to a region of interest. |
| `do_3D` | `true` to enable 3D segmentation. Defaults to `false`. |

### Conversion columns (`pymif 2zarr`)

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
| `shards` | `--shards` |
| `shard_target_mb` | `--shard_target_mb` |
| `shard_exclude_axes` | `--shard_exclude_axes` |

`chunk_size` takes precedence over `max_size(MB)` when both are set. `shards`
only applies with `zarr_format` `3` (`auto` sizes shards automatically,
default target `64` MB via `shard_target_mb`; `shard_exclude_axes` defaults to
`t c`, use `none` to allow every axis to merge). Paths containing commas must
be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,cellpose_diameter,cellpose_channels,do_3D,labels_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 1 1 512 512,2,30,0,False,/scratch/user/labels/sample_01_labels.ome.zarr
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,,,2,25,0,False,
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--default_outdir` | `results` | Base directory for `output`/`labels_output` when a row leaves them blank. |
| `--resolution_level` | `0` | Pyramid level when a row has none. |
| `--cellpose_halo` | `20` | Segmentation-context halo when a row has none. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_cellpose/main.nf \
    -lib /path/to/mif-resources/lib \
    --input_csv /path/to/input.csv \
    -resume
```

For batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit:

```bash
mkdir -p logs
sbatch template_run.sh
```
