# complete_conv_tiled_cellpose

`complete_conv_tiled_cellpose` converts microscopy datasets to OME-Zarr, splits
each converted image into overlapping tiles, segments every tile with
Cellpose-SAM, and stitches the per-tile segmentations into a single, globally
consistent label image (OME-Zarr).

It is the "complete" version of `conv_tiled_cellpose`: that workflow stops after
segmenting the tiles, whereas this one also matches labels across tile borders
and writes the merged mask.

## Steps

1. **`PYMIF_CONVERSION`** – converts each input dataset to OME-Zarr with
   `pymif 2zarr`. If the target store already exists, conversion is skipped.
2. **`TILE_CONSTRUCTOR`** – computes the tile grid (and optional ROI) for each
   converted image and writes a tile CSV.
3. **`PREALLOCATE_ZARR`** – pre-allocates the merged `uint32` label store.
4. **`CELLPOSE_SAM_ZARR_MIF`** – segments every tile on a GPU. Tiles whose
   segmentation was already published (see [Tile scratch](#tile-scratch)) are
   skipped.
5. **`FIND_TILE_NEIGHBORS`** / **`LOCAL_LABEL_MATCHER`** – find the right,
   bottom and behind neighbour of each tile and match overlapping labels between
   neighbouring tiles.
6. **`LABEL_SOLVER`** – solves the label-matching graph of each dataset
   (multicut) into a global label mapping.
7. **`WRITE_RELABELED_TILE`** – writes every tile into the merged label store
   using the global mapping.
8. **`CLEANUP_TILE_SCRATCH`** – removes the per-tile segmentations once the
   merged mask is complete.

## Input CSV

The CSV must contain one dataset per row. The only required column is:

- `input`: path to the input microscopy dataset.

Empty cells fall back to the defaults listed below. Paths containing commas
must be CSV-quoted.

### General columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Path of the converted OME-Zarr. Defaults to `<outdir>/zarr/<dataset_id>.ome.zarr`. Use an absolute path. |
| `labels_output` | Path of the merged label OME-Zarr. Defaults to `<outdir>/masks/<dataset_id>_merged_mask.ome.zarr`. Use an absolute path. |

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

When both `chunk_size` and `max_size(MB)` are present, `chunk_size` takes
precedence. `shards` is only valid with `zarr_format` `3`: use `auto` to size
shards automatically, or an explicit shard shape matching the dataset axes
(e.g. `1 1 8 2160 4096`). `shard_target_mb` only applies when `shards` is `auto`
(default `64`). `shard_exclude_axes` defaults to `t c`; use `none` to allow every
axis to merge.

### Tiling columns

| CSV column | Purpose |
| --- | --- |
| `tile_size_x`, `tile_size_y`, `tile_size_z` | Tile size along each axis. If omitted, derived from the chunk shape of the store at `resolution_level`. |
| `tile_overlap` | Overlap between neighbouring tiles (same units as the tile size). Defaults to `0`. |
| `resolution_level` | Pyramid level used for tiling and segmentation. Defaults to `--resolution_level` (0). |
| `use_physical_units` | `true` to interpret tile sizes, overlap and ROI bounds in physical units instead of pixels. Defaults to `false`. |
| `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` | Restrict processing to a region of interest. Defaults to the full image. |

### Cellpose columns

| CSV column | Purpose |
| --- | --- |
| `cellpose_diameter` | Expected object diameter in pixels (integer). |
| `cellpose_niter` | Number of Cellpose iterations (integer). |
| `cellpose_channels` | Channels used for segmentation. |
| `cellpose_flow_threshold` | Cellpose flow threshold (float). |
| `cellpose_cellprob_threshold` | Cellpose cell-probability threshold (float). |
| `cellpose_halo` | Extra segmentation context around each tile, in the same units as `tile_overlap`. Defaults to `20`. |
| `do_3D` | `true` to enable 3D segmentation. Defaults to `false`. |

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,num_levels,tile_size_x,tile_size_y,tile_overlap,resolution_level,cellpose_diameter,cellpose_channels,cellpose_halo,do_3D,labels_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 1 1 512 512,2,4,2048,2048,128,0,30,0,20,False,/scratch/user/masks/sample_01_mask.ome.zarr
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,,,2,4,,,,0,25,0,,False,
```

## Parameters

These apply to the whole run, not to a single image.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--outdir` | `results` | Base directory for the default converted-zarr and mask locations. |
| `--resolution_level` | `0` | Pyramid level used when a row has no `resolution_level`. |
| `--label_offset_step` | `100000` | Label-ID offset between tiles (`tile_index * label_offset_step`). Must exceed the largest number of labels in any single tile, otherwise IDs collide. |
| `--multicut_beta` | `0.5` | `LABEL_SOLVER` trade-off for cutting edges of the label-matching graph. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. If unset, the default Cellpose-SAM model is used. |

## Outputs

```text
<output column>                                      converted OME-Zarr  (default: <outdir>/zarr/<dataset_id>.ome.zarr)
<labels_output column>                               merged label mask   (default: <outdir>/masks/<dataset_id>_merged_mask.ome.zarr)
```

### Tile scratch

Per-tile segmentations are published (hard-linked) to
`/scratch/$USER/cellpose_tiles_published/<dataset_id>/`. On `-resume`, or on a
re-run, tiles that already have a non-empty store there are not segmented again.
Once the merged mask of a dataset is complete, this directory is deleted.

Because the location is fixed to `/scratch/$USER`, `$USER` must be set and
`/scratch` must be writable. Two runs of the same `dataset_id` at the same time
share (and can corrupt) this directory.

## Running the workflow

The supplied configuration targets the EMBL Heidelberg Slurm cluster, uses
Apptainer for its containers and requests GPUs for the Cellpose step.

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/complete_conv_tiled_cellpose/main.nf \
    -lib /path/to/mif-resources/lib \
    --input_csv /path/to/input.csv \
    --outdir /scratch/$USER/complete_conv_tiled_cellpose/results \
    -resume
```

For a batch submission, update the marked values in `template_run.sh`, create
the launch and log directories, and submit the script:

```bash
mkdir -p /scratch/$USER/complete_conv_tiled_cellpose/logs
cd /scratch/$USER/complete_conv_tiled_cellpose
sbatch /path/to/mif-resources/workflows/complete_conv_tiled_cellpose/template_run.sh
```

`sbatch` writes the Slurm logs to `logs/` relative to the directory it is
submitted from, so submit from the same directory that contains `logs/`.
