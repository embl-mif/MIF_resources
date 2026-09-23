# conv_mip_label_filter_mask_cellpose

Runs each dataset through seven stages:

1. Convert to OME-Zarr with `pymif 2zarr`.
2. Build a multiscale 2D maximum-intensity projection (MIP) over Z.
3. Otsu-threshold one MIP channel (`otsu_channel`) and label it with
   `skimage.measure.label`, rebuilding the full multiscale pyramid (same
   number of levels, chunksize and downsample factors as the MIP) so label
   IDs stay consistent across every resolution.
4. Measure every Otsu-labeled object against the MIP (every channel by
   default) with `skimage.measure.regionprops`.
5. Keep only the objects whose measured row satisfies the requested area
   (and/or any other regionprops-based) limits, and relabel the survivors to
   a contiguous `1..N` id space.
6. Zero `clean_channel` of the MIP outside the filtered objects, and append
   it as a new channel of a new OME-Zarr (every original MIP channel plus
   this one).
7. Segment that new OME-Zarr with Cellpose-SAM, using the appended clean
   channel together with any other channel as `cellpose_channels`.

## Input CSV

One dataset per row. Required columns:

- `input`: path to the input microscopy dataset.
- `otsu_channel`: zero-based index, in the MIP output, of the channel to
  Otsu-threshold and label.

### Output columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names and default output names. Defaults to `<input basename>_<row number>`. |
| `output` | Name of the converted OME-Zarr. Must be absolute if set. Defaults to `<default_outdir>/zarr/<dataset_id>.ome.zarr`. |
| `mip_output` | Name of the MIP OME-Zarr, same rule. Defaults to `<default_outdir>/mip/<dataset_id>_mip.ome.zarr`. |
| `otsu_output` | Name of the Otsu label OME-Zarr, same rule. Defaults to `<default_outdir>/otsu_labels/<dataset_id>_otsu_labels.ome.zarr`. |
| `features_output` | Name of the (pre-filter) measurements CSV, same rule. Defaults to `<default_outdir>/measurements/<dataset_id>_features.csv`. |
| `filter_output` | Name of the filtered label OME-Zarr, same rule. Defaults to `<default_outdir>/filtered/<dataset_id>_filtered_labels.ome.zarr`. |
| `filter_csv_output` | Name of the filtered (post-filter, relabeled) measurements CSV, same rule. Defaults to `<default_outdir>/filtered/<dataset_id>_filtered_labels.csv`. |
| `clean_output` | Name of the OME-Zarr passed to Cellpose, same rule. Defaults to `<default_outdir>/clean/<dataset_id>_clean.ome.zarr`. |
| `labels_output` | Name of the Cellpose label OME-Zarr, same rule. Defaults to `<default_outdir>/labels/<dataset_id>_labels.ome.zarr`. |

### Otsu-label columns

| CSV column | Purpose |
| --- | --- |
| `otsu_channel` | **Required.** Zero-based channel index (in the MIP output) to Otsu-threshold and label. |
| `otsu_resolution_level` | MIP pyramid level to threshold+label at; every other level is rebuilt from it. Defaults to `resolution_level`, or `--resolution_level` (0) if that is also blank. |
| `otsu_multiplier` | Multiplier applied to the Otsu threshold before binarizing (e.g. `1.2` for a stricter cutoff). Defaults to `1.0`. |

### Measurement columns

| CSV column | Purpose |
| --- | --- |
| `features_resolution_level` | Pyramid level to measure at, for both the MIP and the Otsu labels (they share the same pyramid). Defaults to `otsu_resolution_level`. |
| `features_channels` | MIP channel indices to measure, separated by `;`, `,` or whitespace. Defaults to every channel. |
| `features_extra_properties` | Extra `skimage.measure.regionprops` properties to add to the base structural set, separated by `;`, `,` or whitespace. |
| `features_extra_intensity_properties` | Extra `regionprops` properties that require an intensity image, measured once per channel in `features_channels`. |

Every measured column is suffixed per channel as `<property>_ch<index>` (or
`<property>_ch<index>_<channel name>` when available), e.g.
`mean_intensity_ch1_DAPI`; structural properties are unsuffixed. See
[feature_extraction](../feature_extraction/README.md#measurement-columns) for
the full property reference.

### Filter columns

| CSV column | Purpose |
| --- | --- |
| `filter_min_area` | Minimum object area, in pixels, to keep. Blank means no lower bound. |
| `filter_max_area` | Maximum object area, in pixels, to keep. Blank means no upper bound. |
| `filter_properties` | Extra `property:min:max` bounds on any column of the measurements CSV, separated by `;`. Either bound may be left empty, e.g. `solidity:0.8:;mean_intensity_ch1:100:`. Unknown property names raise a clear error listing what was measured. |

Kept objects are relabeled to a contiguous `1..N` id space; the filtered CSV
keeps the new id as `label` and the pre-filter id as `original_label`.

### Clean-channel column

| CSV column | Purpose |
| --- | --- |
| `clean_channel` | MIP channel index to zero outside the filtered objects. Defaults to `otsu_channel`. |

The clean channel is appended after every original MIP channel (e.g. if the
MIP has channels `0..3`, the clean channel becomes channel `4`). Use that
index together with any other channel in `cellpose_channels` below.

### Cellpose columns

| CSV column | Purpose |
| --- | --- |
| `cellpose_diameter` | Expected object diameter in pixels (integer). |
| `cellpose_niter` | Number of Cellpose iterations (integer). |
| `cellpose_channels` | Channels used for segmentation, e.g. `4;1` to combine the clean channel (index `4`) with channel `1`. |
| `cellpose_halo` | Segmentation-context halo in pixels. Defaults to `20`. |
| `resolution_level` | Pyramid level used for segmentation. Defaults to `--resolution_level` (0). Independent of `otsu_resolution_level`/`features_resolution_level` — the clean-channel OME-Zarr keeps the MIP's full pyramid, so any level can be chosen here. |
| `do_3D` | Any non-empty value enables 3D segmentation; leave blank to disable. |

### Conversion columns (`pymif 2zarr`)

Same set as [conv_cellpose](../conv_cellpose/README.md#conversion-columns-pymif-2zarr):
`microscope`, `chunk_size`, `max_size(MB)`, `scene_index`, `zarr_format`,
`downscale_factor`, `channel_colors`, `channel_names`, `num_levels`, `subset`,
`shards`, `shard_target_mb`, `shard_exclude_axes`. `chunk_size` takes
precedence over `max_size(MB)` when both are set. Paths containing commas must
be CSV-quoted.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,otsu_channel,filter_min_area,filter_max_area,filter_properties,cellpose_diameter,cellpose_channels
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,0,50,5000,form_factor:0.5:,30,4;1
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,0,50,,,25,4
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--default_outdir` | `results` | Base directory for every `*_output` column when a row leaves it blank. |
| `--resolution_level` | `0` | Pyramid level when a row has no `otsu_resolution_level` (or `resolution_level`). |
| `--cellpose_model_path` | none | Custom pretrained Cellpose model. Defaults to Cellpose-SAM. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/conv_mip_label_filter_mask_cellpose/main.nf \
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
