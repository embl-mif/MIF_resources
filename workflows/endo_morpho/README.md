# endo_morpho

Identical to
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md)
through its seven stages (convert, MIP, seed Cellpose-SAM, measure, filter,
clean channel, final Cellpose-SAM), plus two more stages that measure
per-cell endothelial morphology and marker intensities, stripped down from
the standalone `endo_morph` scripts
(`Z:\people\gritti\projects\Candelas_Adrian\endo_morph`):

8. For every final-Cellpose-labelled cell that does not touch the image
   border, compute three nested edge-distance bands measured inward from the
   cell edge (junction, cortical actin, transversal actin), and write them as
   a label OME-Zarr. Cells touching the image border are cut by the field of
   view and are never measured (excluded from the bands entirely).
9. Against those bands, measure per-cell shape morphology (area, perimeter,
   aspect ratio, orientation) and V-cadherin/myosin/actin mean intensity in
   each band (and over the whole cell), and estimate the cell-free ("empty")
   area of the image from the fluorescence signal alone -- never from the
   labels, since Cellpose can both miss cells and over-reach into cell-free
   pixels.

Both new stages always run at the same pyramid level as the final Cellpose
segmentation (`resolution_level`/`meta.resolution_level`), for the same
reason `features_resolution_level` is pinned to `0` for the seed labels in
`conv_mip_cpsam_filt_mask_cpsam`: `CELLPOSE_SAM_ZARR_MIF` only builds pyramid
levels *from* the level it segments onward, so the final labels' own level
`0` is the only one guaranteed to share the MIP's shape.

Unlike the standalone `endo_morph` tool (which analyses one timepoint of one
image/labels pair per invocation), only a single timepoint is supported here
(`endomorph_timepoint`, default `0`) -- the MIPs this workflow produces are
already single-timepoint.

## Input CSV

Every column of
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md#input-csv)
applies unchanged. These columns are added:

### Output columns

| CSV column | Purpose |
| --- | --- |
| `endomorph_regions_output` | Name of the region-bands OME-Zarr (stage 8), same absolute-path-or-blank rule as the other `*_output` columns. Defaults to `<default_outdir>/endomorph_regions/<dataset_id>_endomorph_regions.ome.zarr`. |
| `endomorph_cells_output` | Name of the per-cell CSV (stage 9: morphology + band intensities), same rule. Defaults to `<default_outdir>/endomorph_features/<dataset_id>_endomorph_cells.csv`. |
| `endomorph_summary_output` | Name of the per-image summary CSV (stage 9), same rule. Defaults to `<default_outdir>/endomorph_features/<dataset_id>_endomorph_summary.csv`. |

### Region-bands columns (stage 8, `PYMIF_ENDOMORPH_REGIONS`)

| CSV column | Purpose |
| --- | --- |
| `endomorph_ring_depths` | Three band depths (`R1 R2 R3`), from the cell edge in pixels of the final segmentation's resolution level, separated by `;`, `,` or whitespace. All bands start at the edge, so they are nested (e.g. `6 9 24` -> R1 = 0-6px, R2 = 0-9px, R3 = 0-24px from the edge). Defaults to `6 9 24`. |
| `endomorph_border_margin` | Cells with any pixel within this many pixels of the image border are considered cut by the field of view and are not measured. Defaults to `10`. |

Bands are written to `endomorph_regions_output` as a single-level OME-Zarr
label image: value `0` = not in any band (including cells excluded for
touching the border), `1` = R1, `2` = R2 outside R1, `3` = R3 outside R2.
`PYMIF_ENDOMORPH_REGIONS` writes only this zarr -- no CSV. Downstream,
`PYMIF_ENDOMORPH_FEATURES` treats a cell as measured if it has any painted
band pixel; a cell with none was excluded here for touching the border.

### Feature columns (stage 9, `PYMIF_ENDOMORPH_FEATURES`)

| CSV column | Purpose |
| --- | --- |
| `endomorph_vcad_channel`, `endomorph_myosin_channel`, `endomorph_actin_channel`, `endomorph_nuclei_channel` | MIP channel to use for each marker: either the exact `channel_names` entry (e.g. `Alexa_488`) or its index (e.g. `0`). Left blank, each is matched automatically by keyword against `channel_names` (`488`/`cad` for vcad, `561`/`568`/`myo` for myosin, `647`/`actin` for actin, `hoechst`/`dapi`/`405` for nuclei); an ambiguous or missing match raises a clear error. The `nuclei` channel is only used to exclude it from `endomorph_empty_channels` matching below, not measured itself. |
| `endomorph_empty_channels` | Markers whose signal means cells are present, separated by whitespace, e.g. `myosin actin`. The nuclei stain is excluded by default (it stains only the nucleus, so the rest of the cell is as dark as background in that channel). Defaults to `vcad myosin actin`. |
| `endomorph_empty_smoothing` | Gaussian smoothing of the signal before thresholding, in physical units (e.g. µm). Defaults to `2.0`. |
| `endomorph_empty_min_area` | Empty regions smaller than this (physical units²) are ignored. Defaults to `25.0`. |
| `endomorph_empty_seed_percentile` | An empty region starts where a pixel is darker than this percentile of the cell pixels in *every* `endomorph_empty_channels` channel at once. No such pixel means no empty area (e.g. a confluent field). Defaults to `1.0`. |
| `endomorph_empty_min_contrast` | A channel contributes to the cell signal only if `(cell - background) / cell` reaches this. Defaults to `0.25`. |

See the `endo_morph` README
(`Z:\people\gritti\projects\Candelas_Adrian\endo_morph\README.md`) for the
full algorithm description (band geometry, empty-area estimation) -- the
logic here is ported unchanged, only split across two Nextflow processes and
stripped of the QC (`--save-regions` histogram/plot) outputs.

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,seed_channels,seed_diameter,filter_min_area,filter_max_area,filter_properties,cellpose_diameter,cellpose_channels,endomorph_ring_depths,endomorph_border_margin
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,2,30,50,5000,form_factor:0.5:,30,4;1,6 9 24,10
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,2,30,50,,,25,4,,
```

## Parameters

Same as
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md#parameters)
-- no new workflow-level parameters, every new option is a per-row CSV
column above.

## Outputs

For each dataset (`<id>` = `dataset_id`):

| File | Content |
| --- | --- |
| `endomorph_regions/<id>_endomorph_regions.ome.zarr` | Region-bands label image (stage 8; `0`/`1`/`2`/`3`, see above) |
| `endomorph_features/<id>_endomorph_cells.csv` | One row per measured cell: `label`, centroid, area, perimeter, axis lengths, `aspect_ratio`, `orientation_deg`, `R1_area_px`/`R2_area_px`/`R3_area_px`, and per-band marker means `R1_vcad_mean`, `R2_myosin_mean`, `R2_actin_mean`, `R3_myosin_mean`, `R3_actin_mean`, `cell_vcad_mean`, `cell_myosin_mean`, `cell_actin_mean` |
| `endomorph_features/<id>_endomorph_summary.csv` | One row for the image: cell counts, empty-area statistics, channel assignment and the settings used |

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/endo_morpho/main.nf \
    -lib /path/to/mif-resources/lib \
    --input_csv /path/to/input.csv \
    -resume
```

For batch submission, update the marked values in `template_run.sh`, create
the log directory, and submit:

```bash
mkdir -p logs
sbatch template_run.sh
```
