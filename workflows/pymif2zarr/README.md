# pymif2zarr

Converts microscopy datasets to OME-Zarr using one `pymif 2zarr` task per row
of an input CSV.

## Input CSV

One dataset per row. Required column:

- `input`: path to the input microscopy dataset.

`output` (path of the converted OME-Zarr) is optional: give an absolute path
to control it, or leave it blank to use
`<default_outdir>/zarr/<dataset_id>.ome.zarr`.

The other columns are optional and map to `pymif 2zarr` arguments:

| CSV column | PyMIF argument |
| --- | --- |
| `dataset_id` | used as `<dataset_id>` in the default `output` path |
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
dataset_id,input,microscope,output,chunk_size,scene_index,zarr_format,downscale_factor,channel_names,num_levels
sample_01,/data/experiment/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 1 1 512 512,0,2,2,"DAPI,GFP",4
sample_02,/data/experiment/sample_02.ome.tiff,opera,,,0,3,1 2 2,"DAPI,GFP",
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--default_outdir` | `results` | Base directory for `output` when a row leaves it blank. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/pymif2zarr/main.nf \
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
