# LABEL_SOLVER

Concatenates the local graphs of all tiles of one source image and solves a single multicut problem, producing the global label mapping used to stitch per-tile segmentations.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `local_graphs` | path | All `local_graph.json` files of the image (from `LOCAL_LABEL_MATCHER`) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `global_mapping` | `tuple(meta, json)` | Global label mapping (`*.json`) |
