# complete_conv_tiled_cellpose

Converts microscopy datasets to OME-Zarr, tiles each converted image, segments
every tile with Cellpose-SAM, then matches and merges labels across tile
borders into a single global label OME-Zarr.

It's the "complete" version of `conv_tiled_cellpose`, which stops after
segmenting the tiles.

## Steps

1. `PYMIF_CONVERSION` – convert to OME-Zarr (`pymif convert`). Skipped if the
   store already exists.
2. `TILE_CONSTRUCTOR` – compute the tile grid (and optional ROI).
3. `PREALLOCATE_ZARR` – pre-allocate the merged `uint32` label store.
4. `CELLPOSE_SAM_ZARR_MIF` – segment every tile on a GPU. Tiles already
   published (see [Tile scratch](#tile-scratch)) are skipped.
5. `FIND_TILE_NEIGHBORS` / `LOCAL_LABEL_MATCHER` – match labels between
   neighbouring tiles.
6. `LABEL_SOLVER` – solve the label-matching graph into a global mapping.
7. `WRITE_RELABELED_TILE` – write every tile into the merged label store.
8. `CLEANUP_TILE_SCRATCH` – remove the per-tile segmentations once done.

## Input CSV

One dataset per row. Required column:

- `input`: path to the input microscopy dataset.

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Path of the converted OME-Zarr. Defaults to `<outdir>/zarr/<dataset_id>.ome.zarr`. |
| `labels_output` | Path of the merged label OME-Zarr. Defaults to `<outdir>/masks/<dataset_id>_merged_mask.ome.zarr`. |
| `tile_size_x`, `tile_size_y`, `tile_size_z` | Tile size per axis. If omitted, derived from the store's chunk shape at `resolution_level`. |
| `tile_overlap` | Overlap between neighbouring tiles. Defaults to `0`. |
| `resolution_level` | Pyramid level for tiling/segmentation. Defaults to `--resolution_level` (0). |
| `use_physical_units` | `true` to interpret sizes/overlap/ROI in physical units. Defaults to `false`. |
| `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` | Restrict processing to a region of interest. |
| `cellpose_diameter`, `cellpose_niter`, `cellpose_channels`, `cellpose_flow_threshold`, `cellpose_cellprob_threshold` | Standard Cellpose-SAM options. |
| `cellpose_halo` | Extra segmentation context around each tile. Defaults to `20`. |
| `do_3D` | `true` to enable 3D segmentation. Defaults to `false`. |

Use absolute paths for `output`/`labels_output`. Conversion columns
(`microscope`, `chunk_size`, `max_size(MB)`, `scene_index`, `zarr_format`,
`downscale_factor`, `channel_colors`, `channel_names`, `num_levels`, `subset`,
`shards`, `shard_target_mb`, `shard_exclude_axes`) map to `pymif convert`
arguments exactly as in [conv_cellpose](../conv_cellpose/README.md). Paths
containing commas must be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,tile_size_x,tile_size_y,tile_overlap,cellpose_diameter,cellpose_channels,do_3D,labels_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 1 1 512 512,2048,2048,128,30,0,False,/scratch/user/masks/sample_01_mask.ome.zarr
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,,,,,,25,0,False,
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--outdir` | `results` | Base directory for `output`/`labels_output` when a row leaves them blank. |
| `--resolution_level` | `0` | Pyramid level when a row has none. |
| `--label_offset_step` | `100000` | Label-ID offset between tiles. Must exceed the largest per-tile label count. |
| `--multicut_beta` | `0.5` | `LABEL_SOLVER` edge-cut trade-off. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Tile scratch

Per-tile segmentations are hard-linked to
`/scratch/$USER/cellpose_tiles_published/<dataset_id>/`. On `-resume` (or a
re-run), already-segmented tiles are skipped; the directory is deleted once
the merged mask for that dataset is complete. Two concurrent runs of the same
`dataset_id` share (and can corrupt) this directory.

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/complete_conv_tiled_cellpose/main.nf \
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
