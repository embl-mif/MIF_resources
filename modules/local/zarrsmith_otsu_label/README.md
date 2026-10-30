# PYMIF_OTSU_LABEL

Otsu-thresholds one channel at a chosen resolution level, labels the mask, and rebuilds the full pyramid by nearest-neighbour resizing, so label ids are consistent across levels. Writes a label OME-Zarr.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, channel / resolution level options, output name) |
| `input_zarr` | path | Intensity OME-Zarr (e.g. the MIP) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | Label OME-Zarr |
| `versions` | path | `versions.yml` |
