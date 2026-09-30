# RUN_ULTRACK

Tracks cells over a whole time-lapse with Ultrack (GPU), from either ilastik foreground/contour probabilities or Cellpose/SAM instance labels. Writes the tracked segments into a preallocated OME-Zarr (Zarr v2 only), plus the tracks table and the lineage graph.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `labels_zarr` | val | Instance-label zarr (or empty when using probabilities) |
| `foreground_zarr`, `contours_zarr` | val | ilastik foreground / contour probability zarrs |
| `raw_zarr` | val | Raw image zarr |
| `output_tracked_zarr` | val | Preallocated zarr the tracked labels are written to |
| `config_toml` | path | Ultrack configuration |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `tracked_zarr` | `tuple(meta, output_tracked_zarr)` | Tracked label zarr |
| `tracks` | `tuple(meta, csv)` | Tracks table |
| `graph` | `tuple(meta, graph)` | Lineage graph |
| `versions` | path | `versions.yml` |
