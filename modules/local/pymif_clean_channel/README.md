# PYMIF_CLEAN_CHANNEL

Zeroes one channel of an intensity OME-Zarr outside the labelled objects and appends the result as a new "clean" channel, in a new OME-Zarr holding all original channels plus this one. Keeps the full pyramid (intensity and labels must share it).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `clean_channel`, `clean_output_name`) |
| `intensity_zarr` | path | Intensity OME-Zarr (e.g. the MIP) |
| `labels_zarr` | path | Label OME-Zarr (e.g. filtered seed labels) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | OME-Zarr with the extra clean channel (`meta.clean_output_name`) |
| `versions` | path | `versions.yml` |
