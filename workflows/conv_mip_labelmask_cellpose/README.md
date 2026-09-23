# conv_mip_labelmask_cellpose

Runs each dataset through five stages:

1. Convert to OME-Zarr with `pymif 2zarr`.
2. Build a multiscale 2D maximum-intensity projection (MIP) over Z.
3. Otsu-threshold one MIP channel (`mask_channel`), label it with
   `skimage.measure.label`, measure every object with `regionprops`, and drop
   objects outside the requested area (and any other regionprops-based)
   bounds.
4. Zero that channel outside the kept objects, and write a new
   single-resolution OME-Zarr holding every original MIP channel plus this
   masked channel appended at the end.
5. Segment that new OME-Zarr with Cellpose-SAM, using the appended masked
   channel together with any other channel as `cellpose_channels`.

## Input CSV

One dataset per row. Required columns:

- `input`: path to the input microscopy dataset.
- `mask_channel`: zero-based index, in the MIP output, of the channel to
  threshold and mask.

### Output columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Name of the converted OME-Zarr. An absolute path is kept as-is; anything else is used as a plain filename and stays in the Nextflow work directory. |
| `mip_output` | Name of the MIP OME-Zarr, same rule as `output`. Defaults to `<dataset_id>_mip.ome.zarr`. |
| `mask_output` | Name of the masked OME-Zarr passed to Cellpose, same rule as `output`. Defaults to `<dataset_id>_masked.ome.zarr`. |
| `labels_output` | Name of the Cellpose label OME-Zarr, same rule as `output`. Defaults to `<dataset_id>_labels.ome.zarr`. |

### Label-mask columns

| CSV column | Purpose |
| --- | --- |
| `mask_channel` | **Required.** Zero-based channel index (in the MIP output) to Otsu-threshold, label, filter and mask. |
| `mask_resolution_level` | MIP pyramid level to read `mask_channel` from. Defaults to `resolution_level`, or `--resolution_level` (0) if that is also blank. |
| `mask_otsu_multiplier` | Multiplier applied to the Otsu threshold before binarizing (e.g. `1.2` for a stricter cutoff). Defaults to `1.0`. |
| `mask_min_area` | Minimum object area, in pixels, to keep. Blank means no lower bound. |
| `mask_max_area` | Maximum object area, in pixels, to keep. Blank means no upper bound. |
| `mask_intensity_channels` | Channel indices (in the MIP output), separated by `;`, `,` or whitespace, to also measure per-object fluorescence intensity on, for use in `mask_filters`. `mask_channel` itself is always measured. |
| `mask_filters` | Extra `property:min:max` bounds on any measured property, separated by `;`. Either bound may be left empty. Available properties include every `skimage.measure.regionprops` shape property (`solidity`, `eccentricity`, `extent`, `major_axis_length`, ...), the derived `form_factor` = `4*pi*area/perimeter^2`, and fluorescence intensity: `mean_intensity`/`min_intensity`/`max_intensity` (for `mask_channel`) and `mean_intensity_ch<c>`/`min_intensity_ch<c>`/`max_intensity_ch<c>` for every channel `c` listed in `mask_intensity_channels`. Example: `solidity:0.8:;mean_intensity_ch2:150:` keeps objects with solidity ≥ 0.8 and a mean intensity of at least 150 in channel `2`. Unknown property names raise a clear error listing what was measured. |

The masked channel is appended after every original MIP channel (e.g. if the
MIP has channels `0..3`, the masked channel becomes channel `4`). Use that
index together with any other channel in `cellpose_channels` below. A
`<mask_output>_regionprops.csv` file with the measurements of every kept
object is written next to the masked OME-Zarr.

### Cellpose columns

| CSV column | Purpose |
| --- | --- |
| `cellpose_diameter` | Expected object diameter in pixels (integer). |
| `cellpose_niter` | Number of Cellpose iterations (integer). |
| `cellpose_channels` | Channels used for segmentation, e.g. `4;1` to combine the masked channel (index `4`) with channel `1`. |
| `cellpose_halo` | Segmentation-context halo in pixels. Defaults to `20`. |
| `do_3D` | Any non-empty value enables 3D segmentation; leave blank to disable. |

Cellpose always segments resolution level `0` of the masked OME-Zarr, since
that OME-Zarr only has a single resolution level (it is not re-pyramided).

### Conversion columns (`pymif 2zarr`)

Same set as [conv_cellpose](../conv_cellpose/README.md#conversion-columns-pymif-2zarr):
`microscope`, `chunk_size`, `max_size(MB)`, `scene_index`, `zarr_format`,
`downscale_factor`, `channel_colors`, `channel_names`, `num_levels`, `subset`,
`shards`, `shard_target_mb`, `shard_exclude_axes`. `chunk_size` takes
precedence over `max_size(MB)` when both are set. Paths containing commas must
be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,mask_channel,mask_min_area,mask_max_area,mask_intensity_channels,mask_filters,cellpose_diameter,cellpose_channels
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,0,50,5000,2,form_factor:0.5:;mean_intensity_ch2:150:,30,4;1
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,0,50,,,,,25,4
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--default_outdir` | `results` | Base directory for `output`/`mip_output`/`mask_output`/`labels_output` when a row leaves them blank. |
| `--resolution_level` | `0` | Pyramid level when a row has no `mask_resolution_level` (or `resolution_level`). |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_mip_labelmask_cellpose/main.nf \
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
