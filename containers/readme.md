
Each folder here is one Docker image. The [Build containers](../.github/workflows/build-containers.yml) GitHub Action builds it and pushes it to `ghcr.io/embl-mif/FOLDER-NAME:VERSION`.

To create a Docker container:

1. Create a new folder `CONTAINER-NAME` (use underscores, the folder name becomes the image name)

1. Inside, create the `Dockerfile` (and any files it copies, e.g. `environment.yml`). Install `procps`, which Nextflow needs

1. Create a file called `VERSION` containing the image tag, e.g. `1.0.0`

1. Open a pull request: the Action builds the image to check it works, without publishing it

1. Merge into `main`: the Action builds the image again and pushes it to `ghcr.io/embl-mif/CONTAINER-NAME:VERSION` (and `:latest`)

1. The first time a new image is published, set its package visibility to public in the GitHub package settings (`Packages` → package → `Package settings` → `Change visibility`)

1. Use it in a module with `container "ghcr.io/embl-mif/CONTAINER-NAME:VERSION"`

To update a container, change its `Dockerfile` and bump `VERSION`. The Action refuses to overwrite a tag that is already published, so that pipelines pinned to that tag keep getting the same image. To rebuild a published tag anyway, or to build an image without changing it, run the workflow manually from the `Actions` tab (`Build containers` → `Run workflow`).

To build and test an image locally:

```bash
docker build -t CONTAINER-NAME:test containers/CONTAINER-NAME
docker run --rm -it CONTAINER-NAME:test python
```
