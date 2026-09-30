# CELLPOSE_SAM_ZARR_MIF

Segments an OME-Zarr image with Cellpose-SAM (GPU) and writes the instance labels as a label OME-Zarr. Retries on GPU out-of-memory / init failures (up to 2 times).

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, plus Cellpose options such as `segmentation_channels`, `diameter`, `halo`, `resolution_level`, `labels_output_name`) |
| `local_image` | path | Local OME-Zarr to segment (or a placeholder file when the image is remote) |
| `image_uri` | val | URI of a remote image (used instead of `local_image`) |
| `pretrained_cellpose_model` | path | Optional custom Cellpose model file |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `labels` | `tuple(meta, zarr)` | Label OME-Zarr (`meta.labels_output_name`) |
