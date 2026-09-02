All image analysis projects developed at the MIF as deposited in the [EMBL Gitlab](https://git.embl.org/).

# SSH access

To be able to clone repositories from Gitlab, one first needs to configure and ssh key. Typically, this is done from within the HPC cluster. To know how to access it, please visit [cluster.md](cluster.md) first.

1. From the terminal, `ssh-keygen -t ed25519 -C "embl_cluster"`
Hit enter to avoid having a password tied to it, the password filed will now just be empty
Enter the following and copy the full text shown: (`/home/username/.ssh/id_ed25519`)

1. Type `cat ~/.ssh/id_ed25519.pub`

1. Navigate to [https://git.embl.org/-/user_settings/ssh keys](https://git.embl.org/-/user_settings/ssh_keys), which you can access by clicking “edit profile” and then “ssh keys” tab.

1. Click on “Add new key”, and copy and paste the full text shown by the previous `cat` command in the “Key” field..

This key will only be valid if you are accessing gitlab from the cluster, name it to indicate it’s to connect from embl cluster in the “Title” field (e.g. “Title” = “embl_cluster”).

The key is valid for 1 year by default. You can change this behavior by modifying the “Expiration date” field.

# Clone repositories

Log into the cluster.

In your home directory, or group directory, if it’s something you want the whole group to have access to, clone your git repository of interest with ssh.

For example, you will see “Code” button on project page, which will have the project’s “address”. For example:

`git clone git@git.embl.org:grp-mif/image-analysis/REPOSITORY_NAME.git`

## Tips

1. Before running the pipeline, make sure you have the latest version of the repository by navigating to the repository folder and typing `git pull`

1.

