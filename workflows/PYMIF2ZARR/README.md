# PYMIF2ZARR

`PYMIF2ZARR` converts microscopy datasets to OME-Zarr using one
`pymif 2zarr` task per row of an input CSV file.

## Input CSV

The CSV must contain one dataset per row. The required columns are:

- `input`: path to the input microscopy dataset.
- `output`: name of the output `.zarr` directory.

The other columns are optional and are converted into the corresponding
`pymif 2zarr` arguments:

| CSV column | PyMIF argument |
| --- | --- |
| `input` | `--input_path` |
| `output` | `--zarr_path` |
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

### Example

```csv
input,microscope,output,chunk_size,max_size(MB),scene_index,zarr_format,downscale_factor,channel_colors,channel_names,num_levels,subset
/data/experiment/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,100,0,2,2,,"DAPI,GFP",4,
/data/experiment/sample_02.ome.tiff,opera,sample_02.ome.zarr,,100,0,3,1 2 2,,"DAPI,GFP",,y=100:2000;x=100:2000
```

For the first row, the conversion module constructs a command equivalent to:

```bash
pymif 2zarr \
    --input_path sample_01.ome.tiff \
    --zarr_path sample_01.ome.zarr \
    --microscope opera \
    --chunk_size 1 1 1 512 512 \
    --scene_index 0 \
    --zarr_format 2 \
    --downscale_factor 2 \
    --channel_names DAPI GFP \
    --num_levels 4
```

When both `chunk_size` and `max_size(MB)` are present, `chunk_size` takes
precedence.

## Running the workflow

From any location:

```bash
nextflow run /path/to/workflows/PYMIF2ZARR/main.nf \
    --input_csv /path/to/input.csv \
```

The CSV is split into a channel containing one element per row. Each element
has the following structure:

```groovy
[meta, row, input_dataset]
```

That channel is passed directly to `PYMIF_CONVERSION`.
