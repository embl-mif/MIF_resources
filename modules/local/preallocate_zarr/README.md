# PREALLOCATE_ZARR

Creates an empty OME-Zarr skeleton (shape, axes, scale, chunks, translation) for the merged output of one source image. Geometry is computed from the tile table's ROI bounds and the input image's metadata.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`) |
| `input_zarr` | val | Source OME-Zarr to inherit metadata from |
| `tiles_csv` | path | Tile table with the ROI bounds |
| `input_image_path` | val | Path of the input image (calibration) |
| `output_zarr` | val | Path of the skeleton to create |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, output_zarr)` | Path of the created store |
| `versions` | path | `versions.yml` |
