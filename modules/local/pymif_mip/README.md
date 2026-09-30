# PYMIF_MIP

Creates a multiscale OME-Zarr holding the maximum-intensity projection over Z of the input, at every resolution level.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `mip_output_name`) |
| `input_zarr` | path | Input OME-Zarr with a Z axis |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, zarr)` | MIP OME-Zarr (`meta.mip_output_name`) |
| `versions` | path | `versions.yml` |
