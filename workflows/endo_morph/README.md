# endo_morph

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
   view and are never measured (excluded from the bands entirely). The same
   zarr also holds the cell-free ("empty") area of the image, estimated from
   the MIP fluorescence signal alone -- never from the labels, since Cellpose
   can both miss cells and over-reach into cell-free pixels -- and the
   not-imaged area.
9. Against those bands, measure per-cell shape morphology (area, perimeter,
   aspect ratio, orientation) and V-cadherin/myosin/actin mean intensity in
   each band (and over the whole cell), and add the image-level empty-area
   statistics from stage 8 to the per-image summary.
10. Merge the per-dataset cells and summary CSVs of all rows into one cells
   CSV and one summary CSV, saved next to the launch `.sh` script.

Both new stages always run at the same pyramid level as the final Cellpose
segmentation (`resolution_level`/`meta.resolution_level`), for the same
reason `features_resolution_level` is pinned to `0` for the seed labels in
`conv_mip_cpsam_filt_mask_cpsam`: `CELLPOSE_SAM_ZARR_MIF` only builds pyramid
levels *from* the level it segments onward, so the final labels' own level
`0` is the only one guaranteed to share the MIP's shape.

## Algorithm

Stages 8 and 9 take the MIP (4 channels: V-cadherin, myosin, nuclei, actin)
and the final Cellpose label image (`0` = background, `1..N` = cells), both at
the same pyramid level and with the same Y/X size. The MIP's metadata gives
the pixel size (so areas and lengths come out in physical units, e.g. µm) and
the `channel_names` used to find which channel holds which marker. A Z axis, if
present, is max-projected, and the time point is selected with
`endomorph_timepoint` (default `0`).

### Bands (R1, R2, R3) -- stage 8

For each cell, the distance of every pixel to the edge of *that cell's* mask is
computed (a mask pixel next to the background has distance 1). Every band
starts at the edge and goes inward to its depth, so the bands are nested:

| Band | Default depth | Pixels included | Markers averaged (stage 9) |
| --- | --- | --- | --- |
| **R1** junction | 6 px | 0-6 px from the edge | V-cadherin |
| **R2** cortical actin | 9 px | 0-9 px from the edge (includes R1) | myosin, actin |
| **R3** transversal actin | 24 px | 0-24 px from the edge (includes R1 and R2) | myosin, actin |
| **cell** | whole mask | all pixels of the cell | V-cadherin, myosin, actin |

Change the depths with `endomorph_ring_depths`. A border shared with a
neighbouring cell counts as an edge (distance is computed per cell, so the
neighbour is background for it). If a cell is thinner than a band, the band
simply covers the whole cell.

Depths are in pixels **of the analysed pyramid level**, so their physical width
changes with it: each pyramid level halves the image size and doubles the pixel
size. To keep the same physical bands at the next level, halve the depths (e.g.
`6 9 24` -> `3 4.5 12`).

The bands are written as one label image: each pixel holds the *narrowest* band
that contains it (`1` = R1, `2` = R2 outside R1, `3` = R3 outside R2). Anything
deeper than R3, background pixels, and cells excluded for touching the image
border are `0`.

### Cells at the image border -- stage 8

Cells touching the image border are cut by the field of view: their shape is
incomplete and the image border is not a real cell edge. They are therefore
**never measured**: they are absent from `_endomorph_cells.csv`, from the mean
values of the summary and from the band label image. They still count in
`n_cells` and as covered (non-empty) area in `empty_area_pct`.

Cellpose masks often stop a few pixels before the image edge, so a cell that
is cut doesn't always reach the last pixel. A cell is therefore treated as
touching the border if its bounding box comes within `endomorph_border_margin`
pixels of any image edge (default **10 px**): a cell is measured only if there
are more than that many background pixels between it and every image edge.
Increase the margin if cut cells still show up in the bands image; decrease it
if whole cells near the edge are dropped. `0` excludes only cells that reach
the very first or last pixel. Like the band depths, the margin is in pixels of
the analysed level.

### Morphology -- stage 9

Computed for every measured cell, with scikit-image `regionprops`:

| Column | Definition |
| --- | --- |
| `area_px`, `area` | Cell area in pixels and in physical units² |
| `perimeter` | Perimeter length (scikit-image estimate) in physical units |
| `major_axis_length`, `minor_axis_length` | Axes of the ellipse with the same second moments as the cell |
| `aspect_ratio` | `major_axis_length / minor_axis_length` (1 = round) |
| `orientation_deg` | Angle of the major axis in degrees, in (-90, 90]: 0 = horizontal, 90 = vertical, positive = counter-clockwise as the image is displayed |

Perimeter and axis lengths assume square pixels (a warning is logged if the
pixel size is anisotropic, and the X size is used).

### Band intensities -- stage 9

For every measured cell, the mean raw pixel value (no background subtraction)
of each marker is computed in the regions listed in the band table above:
`R1_vcad_mean`, `R2_myosin_mean`, `R2_actin_mean`, `R3_myosin_mean`,
`R3_actin_mean`, plus the whole-cell means `cell_vcad_mean`,
`cell_myosin_mean`, `cell_actin_mean`. `R1_area_px`, `R2_area_px` and
`R3_area_px` give the number of pixels in each band (a band is `NaN` for a
marker if it has no pixel). A cell counts as measured if stage 8 painted any
band pixel for it.

### Empty area -- stage 8

The empty area is decided by the fluorescence signal alone, never by the
labels. Cellpose is not perfect in either direction: it can miss cells
(unlabelled pixels that are not empty) and its masks can cover cell-free
pixels (labelled pixels that *are* empty). Both are handled correctly by
looking only at the signal.

Empty pixels are dark in **every** marker channel. Cells are dark in *some*
channel almost everywhere -- cytoplasm has little V-cadherin, a nucleus has
little of any marker -- so the steps below combine the channels. Pixels that
are zero in every channel (e.g. padding outside the acquired field) are **not
imaged** and are left out of all area percentages.

1. **Cell level per channel.** Every channel in `endomorph_empty_channels`
   (default `vcad myosin actin`; the **nuclei stain is excluded**, being dark
   over most of the cell) is smoothed (`endomorph_empty_smoothing`, default 2
   µm, ignoring pixels outside the acquired field). Its *cell level* is the
   median over the interior of the labelled cells (pixels farther than twice
   the smoothing width from the edge of the labelled area). This is the only
   use of the labels: setting a scale, never deciding what is empty.
2. **Seeds.** A pixel seeds an empty region when it is darker than
   `endomorph_empty_seed_percentile` (default 1 %) of the cell pixels **in every
   channel at once**. Nuclei and cytoplasm are dark in one or two channels but
   never in all of them, so they are not seeded. Seeds smaller than
   `endomorph_empty_min_area` are dropped. **No seed means no empty area**,
   which is what a confluent field gives -- there is no threshold to tune for
   that case.
3. **Background level per channel.** Measured where it really is background:
   the median inside the seeds. Taking a low percentile of the whole image
   instead would make the background depend on how much of the field happens
   to be empty, and the measured area would change with the pyramid level.
4. **Cell signal.** `clip((raw - background) / (cell - background), 0, 1)` is 0
   on background and 1 on anything cell-like, bright structures included, and
   is then smoothed. Smoothing *after* the cap keeps the edge of a gap in
   place; smoothing raw intensities first would let bright junctions spill into
   the gap and shrink it. The cell signal is the maximum over the channels, so
   cytoplasm that is dim in one channel is rescued by the others. A channel
   with too little contrast (`endomorph_empty_min_contrast`, default 0.25,
   compared to `(cell - background) / cell`) is left out.
5. **Empty** = the connected regions where the cell signal is below 0.5 (the
   midpoint between background and cell level, i.e. the edge of the region)
   that contain a seed. Labels are not consulted, so a cell-free hole inside a
   cell mask counts as empty. Regions smaller than `endomorph_empty_min_area`
   (default 25 µm²) are ignored.

The summary reports these complementary percentages of the imaged area:

| Column | Meaning |
| --- | --- |
| `empty_area_pct` | No cells, from the signal. **This is the empty area** |
| `unlabelled_area_pct` | No Cellpose label (the naive definition of empty; border cells count as labelled) |
| `undetected_cell_area_pct` | Cell signal but no label: cells missed by Cellpose |
| `empty_in_cells_area_pct` | Empty but covered by a cell mask: segmentation over-reaching |

The last two say how well the segmentation matches the signal in each
direction. Both are usually small; a large value means the masks and the
fluorescence disagree.

| Symptom | What to change |
| --- | --- |
| A field with visible gaps reports 0 % empty | Nothing was dark in every channel at once. Raise `endomorph_empty_seed_percentile` (e.g. 2 or 5), or check whether one channel is bright everywhere (background too high, or a channel that should be left out of `endomorph_empty_channels`) |
| Nuclei or dim cytoplasm are marked empty | Lower `endomorph_empty_seed_percentile` (e.g. 0.5, 0.1): seeds must be pixels that are unambiguously background |
| Cytoplasm inside cells is marked empty | Increase `endomorph_empty_smoothing`, or restrict `endomorph_empty_channels` to the markers that fill the cell (`myosin actin`) |
| Thin protrusions are marked empty | Increase `endomorph_empty_smoothing` |
| Small gaps are missed | Lower `endomorph_empty_smoothing` or `endomorph_empty_min_area` |

Check the `channels` column of the summary (the channel name used for each
marker) before trusting the numbers, and `empty_channel_levels` for the
background and cell level of every channel and which ones were used.

## Input CSV

Every column of
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md#input-csv)
applies unchanged. These columns are added:

### Output columns

| CSV column | Purpose |
| --- | --- |
| `endomorph_regions_output` | Absolute path of the region-bands OME-Zarr (stage 8), or blank. Blank defaults to `<dataset_id>_endomorph_regions.ome.zarr` in the task work directory. |
| `endomorph_cells_output` | Absolute path of the per-cell CSV (stage 9: morphology + band intensities), or blank. Blank defaults to `<dataset_id>_endomorph_cells.csv` in the task work directory. |
| `endomorph_summary_output` | Absolute path of the per-image summary CSV (stage 9), or blank. Blank defaults to `<dataset_id>_endomorph_summary.csv` in the task work directory. |

This workflow has no `default_outdir`: every `*_output` column (these three and
all the inherited ones) is either an absolute path, which is written directly
to that location, or blank, which leaves the file in the work directory of the
process that made it (a relative path is rejected). Set absolute paths for the
outputs you want to keep.

### Region-bands columns (stage 8, `PYMIF_ENDOMORPH_REGIONS`)

| CSV column | Purpose |
| --- | --- |
| `endomorph_ring_depths` | Three band depths (`R1 R2 R3`), from the cell edge in pixels of the final segmentation's resolution level, separated by `;`, `,` or whitespace. All bands start at the edge, so they are nested (e.g. `6 9 24` -> R1 = 0-6px, R2 = 0-9px, R3 = 0-24px from the edge). Defaults to `6 9 24`. |
| `endomorph_border_margin` | Cells with any pixel within this many pixels of the image border are considered cut by the field of view and are not measured. Defaults to `10`. |
| `endomorph_timepoint` | Time index to analyse, for both stages 8 and 9. Defaults to `0` (always `0` for MIPs). |

Bands and empty area are written to `endomorph_regions_output` as a
single-level OME-Zarr label image, at the resolution of the MIP: value `0` =
none of the below (including cells excluded for touching the border), `1` =
R1, `2` = R2 outside R1, `3` = R3 outside R2, `4` = empty area, `5` = not
imaged. Bands keep priority: `4` and `5` are only painted on pixels without a
band, so an empty region overlapping a cell mask shows the band there (the
percentages in the summary, e.g. `empty_in_cells_area_pct`, are computed from
the full empty mask). `PYMIF_ENDOMORPH_REGIONS` also hands the image-level
empty-area statistics to `PYMIF_ENDOMORPH_FEATURES` in a JSON file in its work
directory. Downstream, `PYMIF_ENDOMORPH_FEATURES` treats a cell as measured if
it has any painted band pixel; a cell with none was excluded here for touching
the border.

### Channel and empty-area columns (used by stage 8; the channel columns also by stage 9)

| CSV column | Purpose |
| --- | --- |
| `endomorph_vcad_channel`, `endomorph_myosin_channel`, `endomorph_actin_channel`, `endomorph_nuclei_channel` | MIP channel to use for each marker: either the exact `channel_names` entry (e.g. `Alexa_488`) or its index (e.g. `0`). Left blank, each is matched automatically by keyword against `channel_names` (`488`/`cad` for vcad, `561`/`568`/`myo` for myosin, `647`/`actin` for actin, `hoechst`/`dapi`/`405` for nuclei); an ambiguous or missing match raises a clear error. The `nuclei` channel is only used to exclude it from `endomorph_empty_channels` matching below, not measured itself. |
| `endomorph_empty_channels` | Markers whose signal means cells are present, separated by whitespace, e.g. `myosin actin`. The nuclei stain is excluded by default (it stains only the nucleus, so the rest of the cell is as dark as background in that channel). Defaults to `vcad myosin actin`. |
| `endomorph_empty_smoothing` | Gaussian smoothing of the signal before thresholding, in physical units (e.g. µm). Defaults to `2.0`. |
| `endomorph_empty_min_area` | Empty regions smaller than this (physical units²) are ignored. Defaults to `25.0`. |
| `endomorph_empty_seed_percentile` | An empty region starts where a pixel is darker than this percentile of the cell pixels in *every* `endomorph_empty_channels` channel at once. No such pixel means no empty area (e.g. a confluent field). Defaults to `1.0`. |
| `endomorph_empty_min_contrast` | A channel contributes to the cell signal only if `(cell - background) / cell` reaches this. Defaults to `0.25`. |

See [Algorithm](#algorithm) for what each of these settings does. The logic is
ported unchanged from the standalone `endo_morph` scripts
(`Z:\people\gritti\projects\Candelas_Adrian\endo_morph`), only split across two
Nextflow processes and stripped of the QC histogram plot. The empty-area
labels are part of the regions zarr (the equivalent of `--save-regions`).

### Example

```csv
dataset_id,input,microscope,output,chunk_size,zarr_format,seed_channels,seed_diameter,filter_min_area,filter_max_area,filter_properties,cellpose_diameter,cellpose_channels,endomorph_ring_depths,endomorph_border_margin
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,2,30,50,5000,form_factor:0.5:,30,4;1,6 9 24,10
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,2,30,50,,,25,4,,
```

## Parameters

Same as
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md#parameters)
-- plus one new workflow-level parameter, `--merged_output_dir` (absolute path;
defaults to the launch directory), where the merged CSVs are written.
`template_run.sh` sets it to the directory `sbatch` was run from, i.e. where
the script lives. Every other new option is a per-row CSV column above.

## Outputs

For each dataset (`<id>` = `dataset_id`), at the path given in the matching
`*_output` column, or in the work directory if it is blank:

| File | Content |
| --- | --- |
| `<id>_endomorph_regions.ome.zarr` | Region-bands + empty-area label image (stage 8; `0`-`5`, see above) |
| `<id>_endomorph_cells.csv` | One row per measured cell: `dataset_id`, `label`, centroid, area, perimeter, axis lengths, `aspect_ratio`, `orientation_deg`, `R1_area_px`/`R2_area_px`/`R3_area_px`, and per-band marker means `R1_vcad_mean`, `R2_myosin_mean`, `R2_actin_mean`, `R3_myosin_mean`, `R3_actin_mean`, `cell_vcad_mean`, `cell_myosin_mean`, `cell_actin_mean` |
| `<id>_endomorph_summary.csv` | One row for the image: `dataset_id`, `n_cells` (all labelled cells, including those touching the border), `n_cells_analysed` (rows of the cells CSV), `image_area`, `imaged_area`, the empty-area statistics (`empty_area`, `empty_area_pct`, `unlabelled_area_pct`, `undetected_cell_area_pct`, `empty_in_cells_area_pct`, `empty_threshold`, `empty_seed_area_pct`, `empty_background_ratio`, `empty_channels_used`, `empty_channel_levels`), `mean_area`/`mean_aspect_ratio`/`mean_perimeter` over the measured cells, pixel size and unit, the empty-area settings used, and the `channels` assignment (see [Algorithm](#algorithm)) |

Once all datasets are done, two merged files are written to
`--merged_output_dir` (next to the launch `.sh` script):

| File | Content |
| --- | --- |
| `endomorph_cells_merged.csv` | All `<id>_endomorph_cells.csv` rows, told apart by `dataset_id` |
| `endomorph_summary_merged.csv` | All `<id>_endomorph_summary.csv` rows, one per dataset |

They are overwritten on every run (a `-resume` that reuses cached results
rewrites them from the same inputs).

## Running the workflow

```bash
module load Nextflow/25.10.1

nextflow run /path/to/mif-resources/workflows/endo_morph/main.nf \
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
