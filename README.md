# nf-modules

A collection of modules, subworkflows and workflows for working with
**OME-Zarr** images: converting microscopy data to OME-Zarr, and running
segmentation, feature extraction, pixel classification and tracking on it at scale.
They are written in [Nextflow](https://www.nextflow.io/) and run their tools in
containers.

The resources were developed for the needs of users of the
[Mesoscopic Imaging Facility](https://www.embl.org/groups/mesoscopic-imaging-facility/)
(MIF) at [EMBL Barcelona](https://www.embl.org/sites/barcelona/), but they are
general: anyone can use, adapt or reuse the modules and workflows for their own
OME-Zarr data.

What is in the repository:

- **Nextflow workflows and modules** built around OME-Zarr, see [Workflows](#workflows).
- **Containers** used by those workflows.
- **Guides and configuration for EMBL users**, for getting set up on the EMBL
  Heidelberg HPC cluster and on GitHub.

### Who is this for?

| If you are... | What applies to you |
| --- | --- |
| **Anyone** working with OME-Zarr | The modules in [modules/local/](modules/local/), the subworkflows and the workflows in [workflows/](workflows/), and the conventions below. Use them as they are, or take single modules into your own pipelines. You only need to supply your own Nextflow profile for your compute environment (see [Cluster resources](#cluster-resources)). |
| **An EMBL user** | Additionally, the EMBL Heidelberg HPC setup: [cluster.md](cluster.md), [github.md](github.md), and the `embl_hd` configuration in [conf/](conf/), which already sets SLURM, Apptainer and shared caches. |
| **A user of the MIF in Barcelona** | Additionally, the shared facility locations used in the examples and templates (`/g/mif-hd/shared/...`) and the facility contact under [Help](#help). |

## Getting started

For EMBL users: if you are new to the cluster, follow these guides in order:

1. [cluster.md](cluster.md): get SSH access to the EMBL Heidelberg HPC cluster and log in.
2. [github.md](github.md): set up an SSH key for GitHub on the cluster and clone repositories.

Then clone this repository, for example into your home directory on the cluster:

```bash
git clone git@github.com:embl-mif/nf-modules.git ~/projects/nf-modules
```

If you are outside EMBL, clone it with HTTPS instead
(`git clone https://github.com/embl-mif/nf-modules.git`).

Run `git pull` in that folder before every new analysis so you have the latest
fixes.

## Workflows

Every workflow lives in its own folder under [workflows/](workflows/). Most of them
read an **input CSV with one dataset per row** and process all rows in parallel on
the cluster. Each folder has a `README.md` listing every CSV column and parameter,
a `template_run.sh` for batch submission, and a `nextflow.config` with the
**workflow-specific parameters** (for example the resolution level or the
Cellpose settings). Check it to see what you can set, and
override any parameter on the command line (`--<parameter> <value>`).

### Which workflow do I need?

| I want to... | Workflow |
| --- | --- |
| Convert raw microscopy files to OME-Zarr only | [zarrsmith_convert](workflows/zarrsmith_convert/README.md) |
| Convert and segment whole images with Cellpose-SAM | [conv_cellpose](workflows/conv_cellpose/README.md) |
| Convert and segment images too large for one GPU, tile by tile, and merge labels across tiles | [tiled_cellpose_sam_zarr](subworkflows/local/tiled_cellpose_sam_zarr/main.nf) |
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
| `zarrsmith_convert` | convert | OME-Zarr | no |
| `conv_cellpose` | convert → Cellpose-SAM | OME-Zarr, labels | yes |
| `conv_mip_cellpose` | convert → MIP → Cellpose-SAM | OME-Zarr, MIP, labels | yes |
| `conv_mip_cpsam_filt_mask_cpsam` | convert → MIP → seed Cellpose-SAM → measure → filter → clean channel → final Cellpose-SAM | OME-Zarr, MIP, seed and final labels, measurement CSVs | yes (2 runs) |
| `endo_morph` | as above → edge bands and empty area → morphology and band intensities → merge | as above, plus per-cell and per-image CSVs | yes (2 runs) |
| `org_facs` | convert → tiled Cellpose-SAM → label OME-Zarr → features | OME-Zarr, labels, features CSV | yes |
| `feature_extraction` | measure labels against image | features CSV | no |
| `ilastik_simple` | ilastik pixel classification | probability OME-Zarr | no |
| `ilastik_tile` | tile → ilastik pixel classification | probability OME-Zarr | no |
| `ultrack` | Ultrack tracking from probabilities | tracked labels, tracks table, lineage graph | yes |

The conversion step uses [ZarrSmith](https://github.com/embl-mif/zarrsmith) (`zarrsmith convert`)
and supports the microscopes ZarrSmith knows about (e.g. `luxendo`, `opera`,
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

### 2. Submit as a batch job

Nextflow keeps running for the whole analysis, so we recommend submitting it as a
batch job with the template script we provide. The job keeps running after you log
out or lose your connection.

Each workflow folder has a `template_run.sh`. Copy it to a run directory, edit
the lines marked `# <--- MODIFY` (your e-mail, the paths to `main.nf`, `lib/`,
the work directory and the CSV, and the Cellpose model). Make sure the script
loads both `Nextflow` and `Miniforge3` (`module load Miniforge3/24.1.2-0`).
Miniforge3 is needed because some modules rely on conda environments, which
Nextflow creates on the fly (without it those tasks fail with
`conda: command not found`). Then submit:

```bash
mkdir -p logs
sbatch template_run.sh
```

The job runs the Nextflow controller on the `htc-el8` partition. Nextflow then
submits every task as its own SLURM job, so the controller itself needs little
memory. Logs go to `logs/<job-name>_<job-id>_<user>.out` / `.err`.

If you prefer, you can run the `nextflow run` command from the script inside `tmux`
instead. Either way, `-lib` is **required**: it points Nextflow to the shared helpers
in [lib/Utils.groovy](lib/Utils.groovy), which it doesn't find on its own because
the workflows are in subfolders. From the repository root `-lib lib` is enough.

### 3. Resume after a failure

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
- **Outputs are written directly to their final location, not published.** We
  don't use Nextflow's `publishDir`, which copies results out of the work directory
  after a task finishes. Instead each task writes straight to the absolute output
  path, and leaves a symbolic link in its work directory pointing to it. Nextflow
  follows that link to pass the result to the next step, so large OME-Zarr stores
  are never copied. The link points from the work directory to the output, not
  the other way round, so cleaning the work directory does not delete your results.
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

`org_facs` keeps per-tile segmentations in nextflow `workdir` so tiles that already finished are skipped on a re-run.

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
workflow's `nextflow.config`. `base.config` is generic. `embl_hd.config` is specific
to the EMBL Heidelberg cluster: if you run elsewhere, copy it and adapt the executor,
container and cache settings to your environment.

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
nf-modules/
├── README.md              this file
├── LICENSE                MIT license
├── cluster.md             HPC cluster access
├── github.md              GitHub SSH key and cloning
├── workflows/             runnable pipelines: main.nf, nextflow.config, README.md, template_run.sh
├── subworkflows/local/    reusable groups of steps (tiled Cellpose-SAM, ilastik, Ultrack)
├── modules/local/         single steps (one process each), each with a README describing inputs/outputs
├── conf/                  shared resource and cluster configuration
├── lib/Utils.groovy       CSV parsing and path helpers shared by all workflows (passed with -lib)
├── assets/                empty placeholder files for optional inputs, don't delete
└── containers/            Docker images built for the workflows
```

### Containers

Containers are published to the GitHub Container Registry (`ghcr.io/embl-mif/`).
The recipes in [containers/](containers/) (`zarr_tools`, `cellpose_zarr`,
`multicut_stitching`, `pymif_cpsam`) are built automatically by a GitHub Action.
The `zarr_tools`, `cellpose_zarr` and `multicut_stitching` recipes are adapted from
those of the [EMBL Centre for Bioimage Analysis](https://git.embl.org/grp-cba/containers).
The `zarrsmith` image is built in the [zarrsmith repository](https://github.com/embl-mif/zarrsmith).
The ilastik and Ultrack modules use public images from Docker Hub.

See [containers/readme.md](containers/readme.md) to add or update a container.

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

For questions or bug reports, open an
[issue](https://github.com/embl-mif/nf-modules/issues) on GitHub, or contact the current maintainers.

## Maintainers

This repository is maintained by Nora Vivanco, Data Scientist for Imaging and Omics at the EMBL Data Science Centre, and Nicola Gritti, Bioimaging Project Manager at the Mesoscopic Imaging Facility (MIF), EMBL Barcelona.

If you use any of these resources, please acknowledge the maintainers and their departments.
## License

Copyright (c) 2026 EMBL. Released under the MIT License, see [LICENSE](LICENSE).
