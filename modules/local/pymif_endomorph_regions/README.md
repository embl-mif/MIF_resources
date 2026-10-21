# PYMIF_ENDOMORPH_REGIONS

Builds one label OME-Zarr with, for every cell not touching the image border, N nested edge-distance bands (`1`..`N` = R1..RN, any number N >= 1 of ring depths), plus the empty area (`N+1`) and not-imaged area (`N+2`) estimated from the fluorescence of the given channels. Empty and not-imaged pixels overwrite the bands. Makes no assumption about what the channels stain.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `endomorph_ring_depths`, `endomorph_border_margin`, `endomorph_empty_channels` (default: all channels), `endomorph_empty_*`, resolution levels, output name) |
| `labels_zarr` | path | Final cell labels |
| `intensity_zarr` | path | Intensity OME-Zarr (the MIP, not the cleaned image) |
| `nucleus_zarr` | path | Clean-channel OME-Zarr: its nuclear channel (`meta.endomorph_empty_nucleus_channel`, default the `..._clean` Hoechst/DAPI channel) marks host nuclei, which are never empty |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, regions_zarr)` | Regions label zarr (`meta.endomorph_regions_output_name`) |
| `qc` | `tuple(meta, json)` | QC diagnostics of the empty-area estimation (work dir, not used downstream) |
| `versions` | path | `versions.yml` |
