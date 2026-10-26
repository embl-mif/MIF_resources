# MIF_resources

__[[TOC]]__

This repository collects resources for users of the
[Mesoscopic Imaging Facility](https://www.embl.org/groups/mesoscopic-imaging-facility/)
(MIF) at [EMBL Barcelona](https://www.embl.org/sites/barcelona/):

- **Guides** for getting set up on the EMBL Heidelberg HPC cluster and the EMBL GitLab.
- **Nextflow workflows** for converting microscopy data to OME-Zarr and running
  segmentation, feature extraction, pixel classification and tracking on it at scale.
- **Containers** used by those workflows.

## Getting started

New to the cluster? Follow these guides in order:

1. [cluster.md](cluster.md): get SSH access to the EMBL Heidelberg HPC cluster and log in.
2. [gitlab.md](gitlab.md): set up an SSH key for GitLab on the cluster and clone repositories.

Then clone this repository on the cluster, for example into your home directory:

```bash
git clone git@git.embl.org:grp-mif/image-analysis/mif_resources.git ~/projects/mif_resources
```

Run `git pull` in that folder before every new analysis so you have the latest
fixes.

## Workflows

Every workflow lives in its own folder under [workflows/](workflows/). Most of them
read an **input CSV with one dataset per row** and process all rows in parallel on
the cluster. Each folder has a `README.md` listing every CSV column and parameter,
and a `template_run.sh` for batch submission.

### Which workflow do I need?

| I want to... | Workflow |
| --- | --- |
| Convert raw microscopy files to OME-Zarr only | [pymif2zarr](workflows/pymif2zarr/README.md) |
| Convert and segment whole images with Cellpose-SAM | [conv_cellpose](workflows/conv_cellpose/README.md) |
| Convert and segment images too large for one GPU, tile by tile, and merge labels across tiles | [conv_tiled_cellpose](workflows/conv_tiled_cellpose/README.md) |
| Convert, make a 2D maximum-intensity projection (MIP) and segment the MIP | [conv_mip_cellpose](workflows/conv_mip_cellpose/README.md) |
| Segment a MIP in two passes: seed segmentation, filter objects by size/shape, then a final segmentation guided by the filtered seeds | [conv_mip_cpsam_filt_mask_cpsam](workflows/conv_mip_cpsam_filt_mask_cpsam/README.md) |
| Do the above and also measure endothelial cell morphology, edge-band intensities and cell-free area | [endo_morph](workflows/endo_morph/README.md) |
| Convert, segment tile by tile, and extract per-object features (organoid FACS-style analysis) | [org_facs](workflows/org_facs/README.md) |
| Measure features of existing labels against an existing image | [feature_extraction](workflows/feature_extraction/README.md) |
| Run an ilastik pixel classifier on one OME-Zarr | `ilastik_simple` (see [below](#workflows-without-a-readme)) |
| Run an ilastik pixel classifier tile by tile on several OME-Zarrs | `ilastik_tile` (see [below](#workflows-without-a-readme)) |
| Track cells over time with Ultrack | `ultrack` (see [below](#workflows-without-a-readme)) |

### Overview

| Workflow | Steps | Main outputs | GPU |
| --- | --- | --- | --- |
| `pymif2zarr` | convert | OME-Zarr | no |
| `conv_cellpose` | convert → Cellpose-SAM | OME-Zarr, labels | yes |
| `conv_tiled_cellpose` | convert → tile → Cellpose-SAM per tile → match/merge labels | OME-Zarr, merged labels | yes |
| `conv_mip_cellpose` | convert → MIP → Cellpose-SAM | OME-Zarr, MIP, labels | yes |
| `conv_mip_cpsam_filt_mask_cpsam` | convert → MIP → seed Cellpose-SAM → measure → filter → clean channel → final Cellpose-SAM | OME-Zarr, MIP, seed and final labels, measurement CSVs | yes (2 runs) |
| `endo_morph` | as above → edge bands and empty area → morphology and band intensities → merge | as above, plus per-cell and per-image CSVs | yes (2 runs) |
| `org_facs` | convert → tiled Cellpose-SAM → label OME-Zarr → features | OME-Zarr, labels, features CSV | yes |
| `feature_extraction` | measure labels against image | features CSV | no |
| `ilastik_simple` | ilastik pixel classification | probability OME-Zarr | no |
| `ilastik_tile` | tile → ilastik pixel classification | probability OME-Zarr | no |
| `ultrack` | Ultrack tracking from probabilities | tracked labels, tracks table, lineage graph | yes |

The conversion step uses [PyMIF](https://github.com/grinic/pymif) (`pymif 2zarr`)
and supports the microscopes PyMIF knows about (e.g. `luxendo`, `opera`,
`viventis`, `zeiss`, `scape`, and existing `zarr` stores). It is skipped when
the target OME-Zarr already exists.

## Running a workflow

### 1. Prepare the input CSV

Create a CSV with one row per dataset. Only a handful of columns are required
(usually just `input`); the rest are optional and fall back to defaults when
left blank. A minimal CSV for `conv_cellpose`:

```csv
dataset_id,input,cellpose_diameter,cellpose_channels
sample_01,/g/mif-hd/shared/data/sample_01.ome.tiff,30,0
sample_02,/g/mif-hd/shared/data/sample_02.ome.tiff,25,0
```

The column names differ between workflows. Check the workflow's README for the
full list.

### 2. Run interactively (small tests)

From a login node, with the repository cloned to `~/projects/mif_resources`:

```bash
module load Nextflow/25.10.1

nextflow run ~/projects/mif_resources/workflows/conv_cellpose/main.nf \
    -lib ~/projects/mif_resources/lib \
    --input_csv /path/to/input.csv \
    --cellpose_model_path /g/mif-hd/shared/cellpose_models/cpsam \
    -resume
```

`-lib` is **required**: it points Nextflow to the shared helpers in
[lib/Utils.groovy](lib/Utils.groovy), which it doesn't find on its own because
the workflows are in subfolders.

Nextflow itself stays running for the whole analysis, so for anything longer
than a quick test, submit it as a batch job instead.

### 3. Submit as a batch job (recommended)

Each workflow folder has a `template_run.sh`. Copy it to a run directory, edit
the lines marked `# <--- MODIFY` (your e-mail, the paths to `main.nf`, `lib/`,
the work directory and the CSV, and the Cellpose model), then submit:

```bash
mkdir -p logs
sbatch template_run.sh
```

The job runs the Nextflow controller on the `htc-el8` partition. Nextflow then
submits every task as its own SLURM job, so the controller itself needs little
memory. Logs go to `logs/<job-name>_<job-id>_<user>.out` / `.err`.

### 4. Resume after a failure

Re-submit the same command. With `-resume` (on by default in the templates),
finished tasks are reused from the work directory and only failed or new tasks
run again. Keep the same `-work-dir` between runs, or the cache is lost.

## Common conventions

These apply to all CSV-driven workflows unless their README says otherwise.

- **`dataset_id`** names the tasks and default output files. When blank it
  defaults to `<input basename>_<row number>`. Keep it unique per CSV.
- **Output columns** (`output`, `labels_output`, `features_output`, ...) must be
  **absolute paths**, or blank. A relative path is rejected with an error,
  because it would end up inside the Nextflow work directory.
- **Blank output columns** use a default location. Depending on the workflow
  that is either under `--default_outdir` / `--outdir` (default `results/` in
  the launch directory) or the task's work directory. The workflow's README
  says which. Anything left in the work directory is lost when it's cleaned,
  so set absolute paths for everything you want to keep.
- **Input paths**: use absolute paths.
- **Lists in a cell** (channels, indices) are separated by `;`, `,` or
  whitespace, e.g. `0;2` or `0 2`. A cell that contains a comma must be quoted
  (`"DAPI,GFP"`).
- **Channel indices** count from `0`.
- **Pyramid levels**: `resolution_level` `0` is full resolution, each higher
  level is downsampled. Pixel-based sizes (diameters, depths, margins) refer to
  the level you process.
- **Region of interest**: `x_min`/`x_max`/`y_min`/`y_max`/`z_min`/`z_max`
  restrict processing to a crop, in pixels (or physical units when
  `use_physical_units` is `true`).

### Cellpose-SAM model

The Cellpose workflows use Cellpose-SAM. Pass the shared model with
`--cellpose_model_path /g/mif-hd/shared/cellpose_models/cpsam` (already set in
the templates), or the path to your own fine-tuned model.

### Tiled workflows

`conv_tiled_cellpose` and `org_facs` keep per-tile segmentations in
`/scratch/$USER/cellpose_tiles_published/<dataset_id>/` so tiles that already
finished are skipped on a re-run. Don't run the same `dataset_id` in two runs at
the same time, as they would share that directory. See each workflow's README for
when the directory is cleaned up.

## Workflows without a README

These workflows are still in development. They take their parameters on the
command line rather than from a full CSV, and some settings are hard-coded in
their `main.nf`.

### ilastik_simple

Runs a trained ilastik pixel classifier (`.ilp`) on one OME-Zarr and writes a
multiscale probability OME-Zarr.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--input_zarr` | – (required) | Input OME-Zarr. |
| `--ilastik_project` | – (required) | Trained ilastik project (`.ilp`). |
| `--outdir` | `results` | Output to `<outdir>/probabilities/<dataset_id>_probabilities.zarr`. |
| `--dataset_id` | `test_dataset` | Name used for the output. |

### ilastik_tile

Same as `ilastik_simple`, but for several datasets listed in a CSV, each split
into tiles that are classified in parallel and written into one probability
OME-Zarr per dataset.

Parameters: `--input_csv`, `--ilastik_project`, `--outdir` (default `results`).

| CSV column | Purpose |
| --- | --- |
| `input_image_path` | **Required.** Input OME-Zarr. |
| `dataset_id` | Defaults to `<input basename>_<row number>`. |
| `pixel_class_output` | Output store. Defaults to `<outdir>/pixel_class/<dataset_id>_pixel_class.zarr`. |
| `tile_size_x`, `tile_size_y`, `tile_size_z`, `tile_overlap`, `halo` | Tiling. |
| `resolution_level`, `use_physical_units` | Pyramid level and units. |
| `x_min` ... `z_max`, `t_min`, `t_max` | Region of interest and time range. |

### ultrack

Tracks cells over time with [Ultrack](https://github.com/royerlab/ultrack) on a
GPU, starting from an ilastik probability OME-Zarr.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--probabilities_zarr` | – (required) | ilastik probabilities (e.g. from `ilastik_simple`). |
| `--input_zarr` | – (required) | Raw OME-Zarr. |
| `--ultrack_config` | none | Ultrack TOML configuration. |
| `--outdir` | `results` | Tracked labels go to `<outdir>/tracking/<dataset_id>_tracked.zarr`; the tracks table and lineage graph to `<outdir>/tracking/<dataset_id>/`. |
| `--dataset_id` | `test_dataset` | Name used for the outputs. |

The channels, time range (`t_min`/`t_max`), object size limits and solver
settings are currently hard-coded in [workflows/ultrack/main.nf](workflows/ultrack/main.nf).
Edit them there before a real run. Ultrack runs on an RTX 3090 node, because
its CUDA 12.4 image can't compile for the newest GPUs.

## Cluster resources

Shared resource settings live in [conf/](conf/) and are included by every
workflow's `nextflow.config`:

| File | Content |
| --- | --- |
| [conf/base.config](conf/base.config) | CPU/memory/time per process label (`process_low`, `process_medium`, `process_high`, `process_high_memory`). Failed tasks with exit codes 1, 137, 140 or 143 (often out of memory or time) are retried with proportionally more resources. |
| [conf/embl_hd.config](conf/embl_hd.config) | EMBL Heidelberg profile: SLURM executor, Apptainer and conda caches under `/g/mif-hd/shared/`, partition choice (`htc-el8`, `bigmem` for high memory per CPU, `gpu-el8` for GPU tasks). |
| [conf/modules.config](conf/modules.config) | Per-process overrides, e.g. the GPU type and publish directory for Ultrack. |

Each workflow's `nextflow.config` caps requests at `max_memory = 512.GB`,
`max_cpus = 48` and `max_time = 14.d`. Override them on the command line
(e.g. `--max_memory 256.GB`) if needed.

Containers are pulled automatically by Apptainer and cached in
`/g/mif-hd/shared/cache_apptainer/`. The first run of a new container can take a
while.

## Repository layout

```text
mif_resources/
├── README.md              this file
├── cluster.md             HPC cluster access
├── gitlab.md              GitLab SSH key and cloning
├── workflows/             runnable pipelines: main.nf, nextflow.config, README.md, template_run.sh
├── subworkflows/local/    reusable groups of steps (tiled Cellpose-SAM, ilastik, Ultrack)
├── modules/local/         single steps (one process each), each with a README describing inputs/outputs
├── conf/                  shared resource and cluster configuration
├── lib/Utils.groovy       CSV parsing and path helpers shared by all workflows (passed with -lib)
├── assets/                empty placeholder files for optional inputs, don't delete
└── containers/            Docker images built for the workflows
```

### Building containers

See [containers/readme.md](containers/readme.md) to build a container and push it
to the EMBL GitLab registry.

### Adding a workflow

1. Create `workflows/<name>/` with `main.nf` and a `nextflow.config` that
   includes the three files in `conf/` (copy an existing one).
2. Reuse modules from `modules/local/` and subworkflows from
   `subworkflows/local/`. Add a new module with its own `README.md` if needed.
3. Parse the CSV with the helpers in `lib/Utils.groovy` (e.g.
   `Utils.resolveOutput` for output columns) so the conventions above hold.
4. Add a `README.md` (CSV columns, parameters, how to run) and a
   `template_run.sh`, and add the workflow to the tables in this README.

## Help

For questions, contact the
[Mesoscopic Imaging Facility](https://www.embl.org/groups/mesoscopic-imaging-facility/).
