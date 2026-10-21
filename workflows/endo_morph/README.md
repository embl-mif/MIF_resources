# endo_morph

Identical to
[conv_mip_cpsam_filt_mask_cpsam](../conv_mip_cpsam_filt_mask_cpsam/README.md)
through its seven stages (convert, MIP, seed Cellpose-SAM, measure, filter,
clean channel, final Cellpose-SAM), plus two more stages that measure
per-cell endothelial morphology and marker intensities, stripped down from
the standalone `endo_morph` scripts
(`Z:\people\gritti\projects\Candelas_Adrian\endo_morph`):

8. For every final-Cellpose-labelled cell that does not touch the image
   border, compute N nested edge-distance bands measured inward from the
   cell edge (any number of `endomorph_ring_depths`), and write them as
   a label OME-Zarr. Cells touching the image border are cut by the field of
   view and are never measured (excluded from the bands entirely). The same
   zarr also holds the cell-free ("empty") area of the image, estimated from
   the fluorescence of `endomorph_empty_channels` alone -- never from the
   labels, since Cellpose can both miss cells and over-reach into cell-free
   pixels -- and the not-imaged area.
9. Against that regions zarr, measure per-cell shape morphology (area,
   perimeter, aspect ratio, orientation) and the mean intensity of every
   measured channel in each band (and over the whole cell), and write the
   per-image summary (cell counts and empty-area statistics).
10. Merge the per-dataset cells and summary CSVs of all rows into one cells
   CSV and one summary CSV, saved next to the launch `.sh` script.

Both new stages always run at the same pyramid level as the final Cellpose
segmentation (`resolution_level`/`meta.resolution_level`), for the same
reason `features_resolution_level` is pinned to `0` for the seed labels in
`conv_mip_cpsam_filt_mask_cpsam`: `CELLPOSE_SAM_ZARR_MIF` only builds pyramid
levels *from* the level it segments onward, so the final labels' own level
`0` is the only one guaranteed to share the MIP's shape.

## Algorithm

Stages 8 and 9 take the MIP (any number of channels) and the final Cellpose
label image (`0` = background, `1..N` = cells), both at the same pyramid level
and with the same Y/X size. The MIP's metadata gives the pixel size (so areas
and lengths come out in physical units, e.g. µm) and the `channel_names` that
`endomorph_empty_channels` and `endomorph_measure_channels` refer to. A Z axis,
if present, is max-projected, and the time point is selected with
`endomorph_timepoint` (default `0`).

### Bands (R1 ... RN) -- stage 8

For each cell, the distance of every pixel to the edge of *that cell's* mask is
computed (a mask pixel next to the background has distance 1). Every band
starts at the edge and goes inward to its depth, so the bands are nested. The
number of bands is the number of `endomorph_ring_depths` (default `6 9 24`,
i.e. three):

| Band | Default depth | Pixels included |
| --- | --- | --- |
| **R1** | 6 px | 0-6 px from the edge |
| **R2** | 9 px | 0-9 px from the edge (includes R1) |
| **R3** | 24 px | 0-24 px from the edge (includes R1 and R2) |
| **cell** | whole mask | all pixels of the cell |

Change the depths with `endomorph_ring_depths`. A border shared with a
neighbouring cell counts as an edge (distance is computed per cell, so the
neighbour is background for it). If a cell is thinner than a band, the band
simply covers the whole cell.

Depths are in pixels **of the analysed pyramid level**, so their physical width
changes with it: each pyramid level halves the image size and doubles the pixel
size. To keep the same physical bands at the next level, halve the depths (e.g.
`6 9 24` -> `3 4.5 12`).

The bands are written as one label image: each pixel holds the *narrowest* band
that contains it (`1` = R1, `2` = R2 outside R1, ... `N` = RN outside R(N-1)).
Anything deeper than RN, background pixels, and cells excluded for touching the
image border are `0`; the empty area is `N+1` and the not-imaged area `N+2`.

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
of each channel in `endomorph_measure_channels` (default: all) is computed in
every band and in the whole cell: `R1_<channel>_mean` ... `RN_<channel>_mean`
and `cell_<channel>_mean`. `R1_area_px` ... `RN_area_px` give the number of
pixels in each band (the mean is `NaN` for a band with no pixel). A cell counts
as measured if stage 8 painted any band pixel for it.

### Empty area -- stage 8

The empty area is decided by the fluorescence signal alone, never by the
labels. Cellpose is not perfect in either direction: it can miss cells
(unlabelled pixels that are not empty) and its masks can cover cell-free
pixels (labelled pixels that *are* empty). Both are handled correctly by
looking only at the signal.

Empty pixels are dark in **every** empty-area channel. Cells are dark in *some*
channel almost everywhere -- e.g. the cytoplasm has little of a junction
marker, a nucleus has little of any cytoplasmic marker -- so the steps below
combine the channels. Pixels that
are zero in every channel (e.g. padding outside the acquired field) are **not
imaged** and are left out of all area percentages.

1. **Cell level per channel.** Every channel in `endomorph_empty_channels`
   (leave out a channel that is dark over most of the cell, such as a nuclei stain) is smoothed (`endomorph_empty_smoothing`, default 2
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
| Cytoplasm inside cells is marked empty | Increase `endomorph_empty_smoothing`, or restrict `endomorph_empty_channels` to the channels that fill the cell |
| Thin protrusions are marked empty | Increase `endomorph_empty_smoothing` |
| Small gaps are missed | Lower `endomorph_empty_smoothing` or `endomorph_empty_min_area` |

Check `<id>_endomorph_empty_qc.json` (in the work directory of
`PYMIF_ENDOMORPH_REGIONS`) before trusting the numbers: `empty_channel_levels`
gives the background and cell level of every channel and which ones were used.

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
| `endomorph_ring_depths` | Band depths (`R1 R2 ... RN`, any number of bands >= 1, strictly increasing), from the cell edge in pixels of the final segmentation's resolution level, separated by `;`, `,` or whitespace. All bands start at the edge, so they are nested (e.g. `6 9 24` -> R1 = 0-6px, R2 = 0-9px, R3 = 0-24px from the edge; `3 8` -> two bands). Defaults to `6 9 24`. |
| `endomorph_border_margin` | Cells with any pixel within this many pixels of the image border are considered cut by the field of view and are not measured. Defaults to `10`. |
| `endomorph_timepoint` | Time index to analyse, for both stages 8 and 9. Defaults to `0` (always `0` for MIPs). |

With `N` bands, bands and empty area are written to `endomorph_regions_output` as a
single-level OME-Zarr label image, at the resolution of the MIP: `1`..`N` =
R1..RN (each pixel holds the narrowest band containing it), `N+1` = empty
area, `N+2` = not imaged, `0` = none of these (including cells excluded for
touching the border). Empty and not-imaged pixels overwrite the bands, so a
cell pixel that is empty is not counted in that cell's bands.
`PYMIF_ENDOMORPH_FEATURES` reads everything it needs, including the empty
area, from this zarr, and treats a cell as measured if it has any band pixel;
a cell with none was excluded here for touching the border.

### Channel and empty-area columns (stage 8, `PYMIF_ENDOMORPH_REGIONS`; the measured channels are stage 9)

Channels are given as the exact `channel_names` entry of the MIP (e.g. `Alexa_488`, no spaces or commas in the name) or its index (e.g. `0`), separated by whitespace, `,` or `;`. Nothing in the workflow assumes what a channel stains.

| CSV column | Purpose |
| --- | --- |
| `endomorph_empty_channels` | Channels whose signal means cells are present. Defaults to all channels; leave out one that stains only part of the cell, such as a nuclei stain (the rest of the cell is as dark as background in that channel, which makes the empty area smaller). |
| `endomorph_measure_channels` | Stage 9: channels whose mean is computed in every band and in the whole cell. Defaults to all channels. |
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
dataset_id,input,microscope,output,chunk_size,zarr_format,seed_channels,seed_diameter,filter_min_area,filter_max_area,filter_properties,cellpose_diameter,cellpose_channels,endomorph_ring_depths,endomorph_border_margin,endomorph_empty_channels,endomorph_measure_channels
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,opera,sample_01.ome.zarr,1 1 1 512 512,2,2,30,50,5000,form_factor:0.5:,30,4;1,6 9 24,10,Alexa_488 Alexa_568 Alexa_647,
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,opera,sample_02.ome.zarr,,2,2,30,50,,,25,4,3 8,,0 1 3,0 1 3
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
| `<id>_endomorph_regions.ome.zarr` | Region-bands + empty-area label image (stage 8; `1`..`N` = bands, `N+1` = empty, `N+2` = not imaged, `0` = the rest) |
| `<id>_endomorph_cells.csv` | One row per measured cell: `dataset_id`, `label`, centroid, area, perimeter, axis lengths, `aspect_ratio`, `orientation_deg`, `R<k>_area_px` for every band, and the mean of every measured channel in every band and in the whole cell: `R<k>_<channel>_mean` and `cell_<channel>_mean` (`<channel>` = the channel name, spaces replaced by `_`) |
| `<id>_endomorph_summary.csv` | One row for the image: `dataset_id`, `n_cells` (all labelled cells, including those touching the border), `n_cells_analysed` (rows of the cells CSV), `image_area`, `imaged_area`, the area statistics computed from the regions zarr (`empty_area`, `empty_area_pct`, `unlabelled_area_pct`, `undetected_cell_area_pct`, `empty_in_cells_area_pct`), `mean_area`/`mean_aspect_ratio`/`mean_perimeter` over the measured cells, pixel size and unit, `n_regions` and `measured_channels` (channel name -> index) |

`PYMIF_ENDOMORPH_REGIONS` also writes `<id>_endomorph_empty_qc.json` to its work directory: the QC diagnostics of the empty-area estimation (`empty_channel_levels` = background and cell level of every channel and whether it was used, `empty_threshold`, `empty_background_ratio`, `empty_seed_area_pct`, `empty_channels_used`, and the settings used). Check it before trusting the empty area of a new dataset.

Once all datasets are done, two merged files are written to
`--merged_output_dir` (next to the launch `.sh` script):

| File | Content |
| --- | --- |
| `endomorph_cells_merged.csv` | All `<id>_endomorph_cells.csv` rows, told apart by `dataset_id` |
| `endomorph_summary_merged.csv` | All `<id>_endomorph_summary.csv` rows, one per dataset |

They are overwritten on every run (a `-resume` that reuses cached results
rewrites them from the same inputs).

## Optional: gaps (cell-free area)

`PYMIF_ENDOMORPH_REGIONS` can also quantify the gaps with a finer method than the empty area of stage 8, into
**separate outputs**. The regions zarr, its empty label (`N+1`), the empty-area QC and the `empty_area_*` summary columns
are computed exactly as without it. It is **off unless `--gap_channels` is given**, and it is controlled by workflow
parameters, not CSV columns, so the CSV (part of every task's input) stays unchanged:

```bash
nextflow run .../workflows/endo_morph/main.nf -lib .../lib --input_csv input.csv \
    --gap_channels "Alexa_488 Alexa_647" \
    --gap_nucleus_channel HOECHST_33342_clean \
    --gaps_outdir /absolute/path/results_gaps \
    -resume
```

- `--gap_channels`: the cell-body markers (e.g. VE-cadherin and actin) of the MIP; a gap must be dark in all of them.
  Not the nuclear stain.
- `--gap_nucleus_channel`: the nuclear channel of the clean-channel zarr, `<channel>_clean` (the nuclear stain with the
  parasite DNA removed, e.g. `HOECHST_33342_clean`; it must match exactly). Its bright pixels are host nuclei and are
  never gap; parasite DNA is not excluded. Default `auto` (the channel named like hoechst/dapi, preferring `..._clean`);
  `none` = no nucleus exclusion.
- `--gaps_outdir`: absolute output folder. Per dataset: `<id>_gaps.ome.zarr` (1 = gap, 2 = not imaged),
  `<id>_gaps_summary.csv` (`gap_area`, `gap_pct`, `n_gaps`, `largest_gap_area`, ...), `<id>_gaps_qc.json`; all summaries
  in `gaps_summary_merged.csv`.
- Optional: `--gap_extra_args "--gap-smoothing 0.25 ..."`, any other `--gap-*` option of `create_endomorph_regions.py`.

The gaps use the same images and resolution levels as the bands (the MIP at `resolution_level`, default 0 = full
resolution, which is what the method was validated on).

**Method.** A gap is dark in every gap channel and is not a host nucleus. Each channel is lightly smoothed (0.3 µm,
`--gap-smoothing`). Its cell level is a low percentile (10, `--gap-cell-percentile`) of the labelled-cell interior, so
ordinary cytoplasm counts as cell; the labels only set this scale. The background is measured in seeds, pixels darker
than the darkest 1 % of the cells in every channel at once. No gap is reported unless that background is clearly darker
than the darkest cells (`--gap-max-floor-ratio`, 0.55), which avoids inventing gaps in a confluent field. A gap is a
region of normalized cell signal below 0.65 (`--gap-edge-threshold`) that contains a core below 0.2
(`--gap-core-threshold`).

**Small regions are not included in the gap percentage.** A gap is counted only if it has a core at background of at
least 1.5 µm² (about 43 px at 0.186 µm/px) and is at least 3 µm² (about 87 px) in total (`--gap-core-min-area`,
`--gap-min-area`). Smaller dark spots (pits of a few pixels, very narrow slits) are not counted, neither in the area nor
in the number of gaps, so `gap_pct` is a slight underestimate (about +0.06 percentage points, +10 % of the gap area, in
one validation image; an upper bound).

Check `<id>_gaps_qc.json` first: `nucleus_channel` should be `HOECHST_33342_clean`, and where a gap exists the
`floor_ratio` of the used channels should be below about 0.55.

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
