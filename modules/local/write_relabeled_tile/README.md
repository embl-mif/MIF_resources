# WRITE_RELABELED_TILE

Relabels one tile's segmentation with the global multicut mapping and writes its non-overlapping core region into the shared merged OME-Zarr store (created by `PREALLOCATE_ZARR`). Runs once per tile, in parallel.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `tile_label` | path | This tile's label zarr |
| `tiles_csv` | path | Tile table |
| `store` | path | Merged output store |
| `global_mapping` | path | Global label mapping (from `LABEL_SOLVER`) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `write_confirmation` | `tuple(meta, done)` | `*.done` marker file |
