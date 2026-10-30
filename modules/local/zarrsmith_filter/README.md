# ZARRSMITH_FILTER

Keeps only the labels whose row in a measurements CSV satisfies the area / property limits, relabels the survivors to `1..N` consistently across all pyramid levels, and writes the filtered label OME-Zarr and the matching filtered CSV.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `filter_min_area`, `filter_max_area`, `filter_properties`, output names) |
| `labels_zarr` | path | Label OME-Zarr to filter |
| `csv_file` | path | Per-label measurements (from `ZARRSMITH_FEATURE_EXTRACTION`) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr_csv` | `tuple(meta, labels_zarr, csv)` | Filtered labels and kept rows (with `original_label`) |
| `versions` | path | `versions.yml` |
