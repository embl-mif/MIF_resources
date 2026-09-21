The High Performance Computing (HPC) cluster is a computational infrastructure operated by EMBL Heidelberg that has several tens of tousands CPUs as well as GPUs available.

The cluster is controlled by a queuing system (SLURM), which manages, prioritizes and runs the submitted jobs depending on resource requestes and pending times.

More information on the HPC can be found [here](https://www.embl.org/internal-information/it-services/hpc-resources-heidelberg/).

# SSH access

To interact with the cluster, one first needs to upload a public ssh key to [pwtools](https://password.embl.org/sshkey). To do so, please follow the instructions at [https://grp-bio-it.embl-community.io/bio-computing-documentation/access-via-ssh/](https://grp-bio-it.embl-community.io/bio-computing-documentation/access-via-ssh/).

# Login

After activating the SSH key, to login from a terminal or cmd, simply type `ssh login1`

Available modules can be discovered and loaded by, e.g., `module avail Nextflow` and `module load Nextflow/24.10.4`.

