# conv_tiled_cellpose_ng

`conv_tiled_cellpose_ng` converts microscopy datasets to OME-Zarr, divides each
dataset into overlapping spatial tiles, and runs Cellpose-SAM independently on
each tile.

The workflow starts one conversion branch per row in an input CSV:

1. `PYMIF_CONVERSION` converts the source dataset with `pymif 2zarr`.
2. `TILE_CONSTRUCTOR` creates a CSV describing the requested tiles.
3. `CELLPOSE_SAM_ZARR_MIF` segments every tile and emits an OME-Zarr label image.

## Input CSV

The CSV must contain one dataset per row.

### Required columns

| Column | Description |
| --- | --- |
| `input` | Path to the input microscopy dataset. |
| `output` | Name of the converted OME-Zarr directory. |
| `tile_size_x` | Tile width along X. |
| `tile_size_y` | Tile height along Y. |

### Dataset and tiling columns

| Column | Required | Description |
| --- | --- | --- |
| `dataset_id` | No | Dataset identifier. Defaults to the input basename. Values should be unique. |
| `tile_size_z` | No | Tile depth along Z. Omit for two-dimensional data. |
| `tile_overlap` | No | Overlap between adjacent tiles. |
| `resolution_level` | No | Zero-based OME-Zarr pyramid level. Defaults to the downstream tool's default when empty. |
| `x_min`, `x_max` | No | Optional X range to tile. |
| `y_min`, `y_max` | No | Optional Y range to tile. |
| `z_min`, `z_max` | No | Optional Z range to tile. |
| `cellpose_channels` | No | Channel indices used for segmentation, for example `0` or `"0 2"`. |
| `labels_output` | No | Requested labels filename. The current `CELLPOSE_SAM_ZARR` module does not yet use this value. |
| `cellpose_diameter` | No | Parsed into metadata, but not forwarded by the current `CELLPOSE_SAM_ZARR` module. |
| `cellpose_niter` | No | Parsed into metadata, but not forwarded by the current `CELLPOSE_SAM_ZARR` module. |
| `t_min`, `t_max` | No | Parsed into metadata, but not used by the current tiling or Cellpose modules. |

Spatial bounds are interpreted as voxel coordinates because the workflow sets
`use_physical_units` to `false`.

The tile halo is configured globally with `--cellpose_halo`; the bundled
configuration sets it to `20`.

### PyMIF conversion columns

These optional columns are passed to `pymif 2zarr`:

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

When both `chunk_size` and `max_size(MB)` are provided, `chunk_size` takes
precedence.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,max_size(MB),scene_index,zarr_format,downscale_factor,channel_colors,channel_names,num_levels,subset,tile_size_x,tile_size_y,tile_size_z,tile_overlap,resolution_level,x_min,x_max,y_min,y_max,z_min,z_max,cellpose_channels,labels_output,cellpose_diameter,cellpose_niter
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,,0,2,2,,"DAPI,GFP",4,,2048,2048,,128,0,,,,,,,"0 2",sample_01_labels.ome.zarr,30,200
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,100,0,2,2,,"DAPI,GFP",4,,2048,2048,,128,0,0,8192,0,8192,,,0,sample_02_labels.ome.zarr,25,150
```

Quote fields containing commas. Empty optional fields are treated as unset.

## Output locations

Converted datasets are written to:

```text
<outdir>/zarr/<output basename>
```

Tile CSV files and Cellpose OME-Zarr outputs remain in the Nextflow work
directory unless `publishDir` rules are added for `TILE_CONSTRUCTOR` and
`CELLPOSE_SAM_ZARR`.

## Running the workflow

The bundled configuration targets the EMBL Heidelberg Slurm cluster and uses
Apptainer for containers.

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_tiled_cellpose_ng/main.nf \
    -work-dir /scratch/$USER/conv_tiled_cellpose_ng/work \
    --input_csv /path/to/batch.csv \
    --outdir /scratch/$USER/conv_tiled_cellpose_ng/results \
    -resume
```

For Slurm submission, copy `template_run.sh`, update the marked values, create
the log and work directories, and submit it:

```bash
mkdir -p logs /scratch/$USER/conv_tiled_cellpose_ng/work
sbatch template_run.sh
```
