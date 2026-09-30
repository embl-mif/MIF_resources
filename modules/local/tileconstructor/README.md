# TILE_CONSTRUCTOR

Splits an OME-Zarr image into a regular grid of (optionally overlapping) tiles and writes their bounds, in voxels as half-open `[min, max)` intervals, to a CSV. Each row has a unique `dataset_id` of the form `<source_dataset_id>__<tile_id>`.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, tile size / overlap options) |
| `image_uri` | val | URI or path of the source OME-Zarr |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `tiles` | `tuple(meta, csv)` | Tile table (`*_tiles.csv`) |
