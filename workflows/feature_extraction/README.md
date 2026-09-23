# feature_extraction

Loads an intensity OME-Zarr and a label OME-Zarr with PyMIF, measures every
label against every requested intensity channel with
`skimage.measure.regionprops`, and writes one CSV of features per dataset.
It is independent of how the labels were produced, so it can be run after
any of this repo's segmentation workflows (`conv_cellpose`,
`conv_mip_cellpose`, `conv_mip_labelmask_cellpose`, `conv_tiled_cellpose`,
`ilastik_*`, `ultrack`, ...), or against any other OME-Zarr label store.

## Input CSV

One dataset per row. Required columns:

- `intensity`: path to the intensity OME-Zarr.
- `labels`: path to the label OME-Zarr (must cover the same field of view as
  `intensity`, and, if present, the same number of timepoints).

### Output columns

| CSV column | Purpose |
| --- | --- |
| `dataset_id` | Identifier used in task names, the `dataset_id` column of the output CSV, and default output names. Defaults to `<intensity basename>_<row number>`. |
| `features_output` | Path of the output CSV. Must be absolute if set. Defaults to `<default_outdir>/measurements/<dataset_id>_features.csv`. |

### Measurement columns

| CSV column | Purpose |
| --- | --- |
| `intensity_resolution_level` | Pyramid level to read from `intensity`. Defaults to `--intensity_resolution_level` (0). |
| `labels_resolution_level` | Pyramid level to read from `labels`. Defaults to `--labels_resolution_level` (0). The spatial (and time) extent at these two levels must match exactly; a mismatch raises a clear error. |
| `channels` | Intensity channel indices to measure, separated by `;`, `,` or whitespace. Defaults to every channel. |
| `extra_properties` | Extra `skimage.measure.regionprops` properties to add to the base structural set (`area`, `bbox`, `centroid`, `perimeter`, `orientation`, `eccentricity`, `solidity`, `extent`, `major_axis_length`, `minor_axis_length`, plus the derived `form_factor` = `4*pi*area/perimeter^2`), measured once per label. Separated by `;`, `,` or whitespace, e.g. `inertia_tensor_eigvals feret_diameter_max`. An unsupported property for the label image's dimensionality raises a clear error. |
| `extra_intensity_properties` | Extra `regionprops` properties that require an intensity image (e.g. `weighted_centroid`), added to the base per-channel set (`mean_intensity`, `min_intensity`, `max_intensity`) and measured once per channel in `channels`. Same separator rules as `extra_properties`. |

Every measured column is suffixed per channel as `<property>_ch<index>` (or
`<property>_ch<index>_<channel name>` when the intensity OME-Zarr has
`channel_names` metadata), e.g. `mean_intensity_ch1_DAPI`. Structural
properties (and `form_factor`) are not channel-suffixed since they only
depend on the label image.

### Example

```csv
dataset_id,intensity,labels,channels,extra_properties,extra_intensity_properties
sample_01,/scratch/user/mip/sample_01_mip.ome.zarr,/scratch/user/labels/sample_01_labels.ome.zarr,0;1,inertia_tensor_eigvals,
sample_02,/scratch/user/mip/sample_02_mip.ome.zarr,/scratch/user/labels/sample_02_labels.ome.zarr,,,
```

## Parameters

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_csv` | – (required) | Input CSV. |
| `--default_outdir` | `results` | Base directory for `features_output` when a row leaves it blank. |
| `--intensity_resolution_level` | `0` | Pyramid level for `intensity` when a row has none. |
| `--labels_resolution_level` | `0` | Pyramid level for `labels` when a row has none. |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/feature_extraction/main.nf \
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
