# FIND_TILE_NEIGHBORS

For every tile of one source image, finds its forward neighbours (right, bottom, behind) from the tile bounds and overlaps, and lists the paths of their segmented label zarrs. Pure metadata operation, no image I/O.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `tiles_csv` | path | Tile table with bounds, overlaps and label-zarr paths |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `tile_neighbors` | `tuple(meta, csv)` | Tile neighbour table (`*.csv`) |
