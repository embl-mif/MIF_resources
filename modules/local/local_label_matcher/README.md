# LOCAL_LABEL_MATCHER

Builds one tile's contribution to the region-adjacency graph used for multicut stitching: boundary edges between its labels and the labels of its right / bottom / behind neighbours (label ids are offset per tile so they are globally unique).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, tile index, overlaps) |
| `self_zarr` | path | This tile's label zarr |
| `right_zarr`, `bottom_zarr`, `behind_zarr` | path | Neighbour label zarrs (placeholders when the neighbour does not exist) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `local_graph` | `tuple(meta, json)` | Local graph (`*.json`), input of `LABEL_SOLVER` |
