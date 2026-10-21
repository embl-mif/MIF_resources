# PYMIF_ENDOMORPH_REGIONS

Builds one label OME-Zarr with, for every cell not touching the image border, N nested edge-distance bands (`1`..`N` = R1..RN, any number N >= 1 of ring depths), plus the empty area (`N+1`) and not-imaged area (`N+2`) estimated from the fluorescence of the given channels. Empty and not-imaged pixels overwrite the bands. Makes no assumption about what the channels stain.

## Input
| Name | Type | Description |
| --- | --- | --- |
| `meta` | val | Sample map (`id`, `endomorph_ring_depths`, `endomorph_border_margin`, `endomorph_empty_channels` (default: all channels), `endomorph_empty_*`, resolution levels, output name) |
| `labels_zarr` | path | Final cell labels |
| `intensity_zarr` | path | Intensity OME-Zarr (the MIP, not the cleaned image) |
| `nucleus_zarr` | path | Clean-channel OME-Zarr; only read for the optional gaps (its `<channel>_clean` nuclear channel excludes host nuclei) |

## Output
| Name | Type | Description |
| --- | --- | --- |
| `zarr` | `tuple(meta, regions_zarr)` | Regions label zarr (`meta.endomorph_regions_output_name`) |
| `qc` | `tuple(meta, json)` | QC diagnostics of the empty-area estimation (work dir, not used downstream) |
| `gaps` | `tuple(meta, gaps_zarr, summary_csv, qc_json)` | Optional, only with `params.gap_channels`: gap label zarr (1 = gap, 2 = not imaged), one-row gap summary and gap QC, copied to `params.gaps_outdir` |
| `versions` | path | `versions.yml` |

## Optional gaps (workflow params)
| Param | Meaning |
| --- | --- |
| `gap_channels` | Cell-body markers, e.g. `"Alexa_488 Alexa_647"`. A gap must be dark in all of them. Gaps are skipped if not set. |
| `gap_nucleus_channel` | Nuclear channel of `nucleus_zarr`, e.g. `HOECHST_33342_clean`. Default `auto`; `none` = no nucleus exclusion. |
| `gaps_outdir` | Absolute output folder, required with `gap_channels`. |
| `gap_extra_args` | Other `--gap-*` options of `create_endomorph_regions.py`, e.g. `"--gap-smoothing 0.25 --gap-min-area 2"`. |

The gaps are separate outputs: the regions zarr and the empty-area QC are the same with or without them.
