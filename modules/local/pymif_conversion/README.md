# PYMIF_CONVERSION

Converts a microscope dataset to OME-Zarr with `pymif 2zarr`. Options (microscope, chunk size, zarr format, subset, channel names/colors, ...) come from the columns of the input CSV row, carried in `meta`.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | CSV row plus `id` and `zarr_output_name` |
| `input_dataset` | path | Input microscopy dataset |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | Converted OME-Zarr (`meta.zarr_output_name`) |
| `versions` | path | `versions.yml` |
