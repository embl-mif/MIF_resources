# PYMIF_FEATURE_EXTRACTION

Measures every label of a label OME-Zarr against an intensity OME-Zarr with `skimage.measure.regionprops` (shape properties plus per-channel intensities) and writes one CSV row per label.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `features_*` options and output name) |
| `intensity_zarr` | path | Intensity OME-Zarr |
| `labels_zarr` | path | Label OME-Zarr |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `csv` | `tuple(meta, csv)` | Feature table (`meta.features_output_name`) |
| `versions` | path | `versions.yml` |
