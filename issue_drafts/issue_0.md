# Issue 0: Installing Snakemake
In order to walk through this learning module, we will create a new Python environment and install Snakemake and other required packages.

Please make sure you have Anaconda installed on your computer before proceeding with the following steps.

If you are working in Windows, you will need to open your Anaconda prompt to follow the steps below. If you are working in Linux or macOS, you can just open your command line.

## Install using environment.yaml
You can create a Conda environment for this tutorial with all the required packages in one step by running the following command:
`
conda env create -f environment.yaml
`

Note: You may have issues installing with the environment.yaml file. This could be due to the fact that Anaconda sometimes has trouble installing Snakemake. If this is the case, follow the instructions below for a manual installation.

## Create environment and install packages manually
You can manually create your Conda environment and install the required packages by following the steps below.

1. Create a new Conda environment: `conda create -n snakemake-tutorial python=3.8`

2. Activate the environment: `conda activate snakemake-tutorial`


3. Install Mamba, a fast and robust replacement for the Conda package manager that is better able to handle the installation of Snakemake: `conda install -c conda-forge mamba`

4. Install Snakemake (please [read the docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more details on installation). If you are working in a Mac or Linux OS, you can install the full version of Snakemake: `mamba install -c conda-forge -c bioconda snakemake`. If you are working in Windows, some of the dependencies in the full version of Snakemake will not work for you. You can install the minimal version of Snakemake:`mamba install -c conda-forge bioconda::snakemake-minimal`

5. Install any additional packages needed for this tutorial (listed in the environment.yaml file) using Conda or Mamba.




**Last Resort**: If neither of the above options work, you can try installing Snakemake or snakemake-minimal with the appropriate Conda command:

`
conda install -c conda-forge bioconda::snakemake
conda install -c conda-forge bioconda::snakemake-minimal
`
