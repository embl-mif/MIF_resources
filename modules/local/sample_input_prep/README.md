# SAMPLE_INPUT_PREP

Scans an experiment directory of OME-TIFF files and generates the sample table (CSV) for the workflows, extracting the microscope model, channel names and pyramid levels from the OME metadata.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `exp_dir` | path | Experiment directory with the images |
| `zarr_ver` | val | Zarr format version to put in the table |
| `sample_table` | val | Name of the CSV to write |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `csv` | path | Generated sample table |
| `versions` | path | `versions.yml` |
