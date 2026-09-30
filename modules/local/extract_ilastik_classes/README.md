# EXTRACT_ILASTIK_CLASSES

Reads an ilastik `.ilp` project and reports its class count and class names as JSON. Run upstream of `PREALLOCATE_ZARR` to set the number of output channels and their labels.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `project` | path | ilastik project file (`.ilp`) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `classes` | path | `ilp_classes.json` with `num_channels` and `channel_labels` |
| `versions` | path | `versions.yml` |
