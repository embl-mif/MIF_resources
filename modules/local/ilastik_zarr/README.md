# ILASTIK_ZARR

Headless ilastik tile worker: reads tiles from an OME-Zarr, runs an ilastik project on them and writes the predictions into an already preallocated OME-Zarr probability store (see `PREALLOCATE_ZARR`).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, optional tile/time ranges `t_min`, `y_min`, `x_min`, ...) |
| `input_zarr` | val | Input OME-Zarr path |
| `master_output_zarr` | val | Preallocated output store to write into |
| `project` | path | ilastik project file (`.ilp`) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, master_output_zarr)` | Path of the store the predictions were written to |
| `versions` | path | `versions.yml` |
