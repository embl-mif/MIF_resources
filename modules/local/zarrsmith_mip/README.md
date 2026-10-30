# ZARRSMITH_MIP

Creates a multiscale OME-Zarr holding the maximum-intensity projection over Z of the input, at every resolution level.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `zarr_output`) |
| `input_zarr` | path | Input OME-Zarr with a Z axis |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | MIP OME-Zarr (`meta.zarr_output`) |
| `versions` | path | `versions.yml` |
