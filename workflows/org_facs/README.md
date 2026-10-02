# org_facs

Converts microscopy datasets to OME-Zarr, tiles each converted image, segments
every tile with Cellpose-SAM, matches and merges labels across tile borders
into an intermediate global mask, converts it to a multiscale label OME-Zarr
with the same metadata as the converted image, and finally measures the labels
against the converted image, writing one features CSV per dataset.

## Steps

1. `PYMIF_CONVERSION` – convert the input to OME-Zarr (`pymif 2zarr`). Skipped
   if the store already exists.
2. `TILED_CELLPOSE_SAM_ZARR` (subworkflow) – segment the converted image tile
   by tile and merge the result into one intermediate label mask, kept only in
   the Nextflow work directory:
   1. `TILE_CONSTRUCTOR` – compute the tile grid (and optional ROI).
   2. `PREALLOCATE_ZARR` – pre-allocate the merged `uint32` label store.
   3. `CELLPOSE_SAM_ZARR_MIF` – segment every tile on a GPU. Tiles already
      published (see [Tile scratch](#tile-scratch)) are skipped.
   4. `FIND_TILE_NEIGHBORS` / `LOCAL_LABEL_MATCHER` – match labels between
      neighbouring tiles.
   5. `LABEL_SOLVER` – solve the label-matching graph into a global mapping.
   6. `WRITE_RELABELED_TILE` – write every tile into the merged label store.
   7. `POPULATE_PYRAMID` – fill pyramid levels 1..N of the merged mask from
      the fully written level 0 (nearest-neighbour `Sample`, so label ids are
      preserved).
3. `PYMIF_CONVERSION` (zarr → zarr) – convert the intermediate mask into the
   final label OME-Zarr (`labels_output_final`) with the metadata of the
   converted OME-Zarr.
4. `PYMIF_FEATURE_EXTRACTION` – measure the label OME-Zarr against the
   converted OME-Zarr (`skimage.measure.regionprops`) and write a features CSV.
   Runs per dataset once its labels are complete.

## Input CSV

One dataset per row. Blank cells fall back to the default listed below. Only
`input` is mandatory. Paths containing commas must be CSV-quoted. A relative
`input` is resolved against the CSV's own directory; every output column
(`output`, `labels_output_final`, `features_output`) must be an **absolute**
path, or blank to use the default.

### General

| CSV column | Required | Default if not set | Description |
| --- | --- | --- | --- |
| `input` | **yes** | – | Path to the input microscopy dataset. |
| `dataset_id` | no | `<input basename>_<row number>` | Identifier used in task names and default output names. Keep it unique: it also names the per-tile scratch directory. |

### 1. Conversion (`pymif 2zarr`)

Convert the input into the OME-Zarr used by every later step. All columns are
optional; unset ones are not passed to `pymif`, so its own defaults apply. The
same `num_levels`, `downscale_factor`, `zarr_format`, ... are reused when the final
label OME-Zarr is built, so labels and image share the pyramid.

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `output` | `<dataset_id>.ome.zarr` in the work dir | Path of the converted OME-Zarr. Absolute path to keep it. Skipped if the store already exists. |
| `microscope` | auto-detected by pymif | Microscope/format of the input (`luxendo`, `opera`, `viventis`, `zeiss`, `zarrv04`, `zarr`, `scape`). |
| `chunk_size` | derived by pymif | Chunk size in `TCZYX` order (or the axes present), space-separated, e.g. `1 1 1 512 512`. Takes precedence over `max_size(MB)`. |
| `max_size(MB)` | pymif default | Maximum chunk size in MB. Ignored if `chunk_size` is set. |
| `scene_index` | pymif default | Scene to convert for multi-scene files (`.czi`). |
| `zarr_format` | pymif default | `2` (NGFF 0.4) or `3` (NGFF 0.5). |
| `num_levels` | pymif default | Number of pyramid levels. |
| `downscale_factor` | pymif default | Pyramid downsampling factor: one value (isotropic) or three (`Z Y X`). |
| `subset` | none (whole image) | Crop before building the pyramid, e.g. `y=10:100:2;x=20:80`. |
| `channel_names` | pymif default | Channel names, space- or comma-separated. |
| `channel_colors` | pymif default | Channel colours (hex or matplotlib names). |
| `shards` | `auto` (zarr v3) | Shard shape, `auto`, or `none` to disable sharding. |
| `shard_target_mb` | pymif default (1024) | Target shard size when `shards` is `auto`. |
| `shard_exclude_axes` | `t c` | Axes never merged into a shard; `none` lets every axis merge. |
| `num_workers` | pymif default | Threads used to write the zarr. |
| `drop_singleton` | `false` (workflow default) | Drop singleton `t`/`c`/`z` axes. Kept `false` here (pymif's own default is `true`) because downstream steps expect the full axis layout. |

See [conv_cellpose](../conv_cellpose/README.md) for more detail on these
columns.

### 2. Tiled Cellpose-SAM segmentation

All columns are optional.

**Tiling and region of interest**

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `tile_size_x`, `tile_size_y`, `tile_size_z` | derived from the store's chunk shape at `resolution_level` | Tile size per axis. |
| `tile_overlap` | `0` | Overlap between neighbouring tiles. |
| `resolution_level` | `--resolution_level` (`0`) | Pyramid level used for tiling and segmentation. |
| `use_physical_units` | `false` | `true` to interpret tile sizes, overlap and ROI in physical units instead of voxels. |
| `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max` | whole image | Restrict processing to a region of interest. Cropped masks cannot be converted by the labels step if the shape no longer matches the image, so avoid combining an ROI with the label conversion unless tested. |

**Cellpose-SAM**

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `cellpose_channels` | all channels | Channel indices to segment (`;`, `,` or space-separated). |
| `cellpose_diameter` | auto-estimated by Cellpose | Expected cell diameter in pixels. |
| `cellpose_niter` | Cellpose default | Number of mask-reconstruction iterations. |
| `cellpose_flow_threshold` | `0.4` | Cellpose flow-error threshold. |
| `cellpose_cellprob_threshold` | `0.0` | Cellpose cell-probability threshold. |
| `cellpose_halo` | `20` | Extra context around each tile (same units as the tile size) so cells at tile borders are segmented properly. |
| `do_3D` | `false` | `true` to force 3D segmentation; otherwise the mode is chosen from the image anisotropy. |

**Final labels**

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `labels_output_final` | `<dataset_id>_merged_mask.ome.zarr` in the work dir | Path of the final label OME-Zarr. Absolute path to keep it. |

The tiles are merged into an *intermediate* label mask (pyramid filled by `POPULATE_PYRAMID`) that is
only kept in the Nextflow work directory. `PYMIF_CONVERSION` (zarr → zarr)
then converts it into `labels_output_final`, a multiscale label OME-Zarr that
reuses the conversion settings (axes, scales, pyramid levels, zarr format) of
the converted image. Do **not** name a CSV column `labels_output`: that would
be picked up as the intermediate mask location.

### 3. Feature extraction

The intensity image is the converted OME-Zarr and the labels are the final
label OME-Zarr; neither needs to be given. All columns are optional.

| CSV column | Default if not set | Description |
| --- | --- | --- |
| `features_output` | `<outdir>/measurements/<dataset_id>_features.csv` | Path of the features CSV. |
| `intensity_resolution_level` | `--intensity_resolution_level` (`0`) | Pyramid level read from the converted image. |
| `labels_resolution_level` | `--labels_resolution_level` (`0`) | Pyramid level read from the labels. The extent at both levels must match. |
| `channels` | every channel | Intensity channel indices to measure (`;`, `,` or space-separated). |
| `extra_properties` | none | Extra `regionprops` properties, e.g. `inertia_tensor_eigvals feret_diameter_max`. |
| `extra_intensity_properties` | none | Extra `regionprops` properties that need an intensity image, e.g. `weighted_centroid`. |

See [feature_extraction](../feature_extraction/README.md) for the full list of
measured properties and the output column naming.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,tile_size_x,tile_size_y,tile_overlap,cellpose_diameter,cellpose_channels,do_3D,labels_output_final,channels,features_output
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,/scratch/user/zarr/sample_01.ome.zarr,1 1 1 512 512,2048,2048,128,30,0,False,/scratch/user/masks/sample_01_mask.ome.zarr,0;1,/scratch/user/features/sample_01.csv
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,,,,,,25,0,False,,,
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--outdir` | `results` | Base directory for the default features CSV (`<outdir>/measurements`). |
| `--resolution_level` | `0` | Pyramid level for tiling/segmentation when a row has none. |
| `--intensity_resolution_level` | `0` | Pyramid level for feature extraction on the image when a row has none. |
| `--labels_resolution_level` | `0` | Pyramid level for feature extraction on the mask when a row has none. |
| `--label_offset_step` | `100000` | Label-ID offset between tiles. Must exceed the largest per-tile label count. |
| `--multicut_beta` | `0.5` | `LABEL_SOLVER` edge-cut trade-off. |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Tile scratch

Per-tile segmentations are hard-linked to
`/scratch/$USER/cellpose_tiles_published/<dataset_id>/`. On `-resume` (or a
re-run), already-segmented tiles are skipped. The directory is not cleaned up
automatically: delete it yourself once the dataset is done. Two concurrent runs of the same
`dataset_id` share (and can corrupt) this directory.

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/org_facs/main.nf \
    -lib /path/to/mif-resources/lib \
    --input_csv /path/to/input.csv \
    --cellpose_model_path /g/mif-hd/shared/cellpose_models/cpsam \
    -resume
```

For batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit:

```bash
mkdir -p logs
sbatch template_run.sh
```
