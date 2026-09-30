# POPULATE_PYRAMID

Fills pyramid levels 1..N of an OME-Zarr from its fully written level 0 (`Average` for continuous data such as probabilities, `Sample` for labels).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `output_zarr` | val | OME-Zarr store to populate (level 0 must already be written) |
| `downsample_method` | val | `Average` or `Sample` |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, output_zarr)` | Path of the populated store |
| `versions` | path | `versions.yml` |
