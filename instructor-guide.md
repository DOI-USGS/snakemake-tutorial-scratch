# How to be an instructor for the Snakemake introduction tutorial

## Prepare
1. Install the conda environment
2. Install mamba if possible

## Get started
- The learner should add you to their repo
- Clone their repo locally
- Checkout each issue branch and test
- Please edit this document as you see fit!

## Issue 0
- Support with installation as needed
- No reviews required

## Issue 1
- Make sure that you've been added as a reviewer - if not, add yourself if possible, or ask the learner to add you.
- Double check that the PR is for merging to their fork, not to USGS-R/snakemake-introduction
- Should be one new file: Snakefile, with these contents
```
rule get_sb_data:
    output:
        "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```
- Make sure that the zip file itself was not committed, just the Snakefile
- Checkout the issue-1 branch locally
- Run `snakemake --cores 1 1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`
- Check that `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip` was created

## Issue 2
- Two changed files: Snakefile and sb_get.py
- New Snakefile contents:
```
rule get_sb_data:
    params:
        sb_item = "5e5d0bb9e4b01d50924f2b36"
    output:
        sb_file = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```
- sb_get.py changes at bottom of file:
```
if __name__ == "__main__":
    sb_item = snakemake.params['sb_item']
    sb_file = snakemake.output['sb_file']
    main(sb_item, sb_file) 
```
- Check out issue-2 branch
- Delete zip file from 1_fetch/tmp
- Run `snakemake --cores 1 1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`
- Check that `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip` was created

## Issue 3
- Change Snakefile and unzip_file.py
- New Snakefile contents at bottom:
```
rule unzip_sb_data:
     input:
         zip_file_path = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
     output:
         "1_fetch/out/tmp/pgdl_nhdhr_{FC091A8F-FC45-46C0-91F9-18379CF0EAAE}_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_69545713_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_80006805_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_86444267_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_86445115_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_105954753_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_107071276_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_107071492_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_107072210_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_111726865_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120018402_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120018788_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120018790_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020150_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020163_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020166_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020167_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020444_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020465_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020466_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020478_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020480_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020497_temperatures.csv",
         "1_fetch/out/tmp/pgdl_nhdhr_120020979_temperatures.csv"
     params: 
         out_dir = '1_fetch/out/'
     script: 
         '1_fetch/unzip_file.py'
```
- Changes to unzip_file.py at bottom:
```
if __name__ == '__main__':
     zip_file_path = snakemake.input['zip_file_path']
     out_dir = snakemake.params['out_dir']
     main(zip_file_path, out_dir)
```
- Check out issue-3 branch
- Run `snakemake --cores 1 unzip_sb_data --delete-all-output`
- Run `snakemake --cores 1 unzip_sb_data`
- Check that unzipped files were created
- Run `snakemake --cores 1 unzip_sb_data --delete-all-output` if you want

## Issue 4
- Change Snakefile and calc_doy_means.py
- New Snakefile contents:
    * New rule all at the top of the file
```
rule all: 
    input: 
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv" 
```
    * New rule: calc_doy_means
```
rule calc_doy_means:
    input:
        in_file = "1_fetch/out/tmp/pgdl_nhdhr_{lake_id}_temperatures.csv"
    output:
        out_file = "2_process/out/doy_{lake_id}.csv"
    script:
        "2_process/calc_doy_means.py"
```
- Changes to calc_doy_means.py at bottom:
```
if __name__ == '__main__': 
    out_file = snakemake.output['out_file']
    in_file = snakemake.input['in_file']
    lake_id = snakemake.wildcards['lake_id']
    main(out_file, in_file, lake_id)
```
- Check out issue-4 branch
- Run `snakemake --cores 1 --delete-all-output` 
    * rule all should be present now, so no need to specify outputs/rules
- Run `snakemake --cores 1`
- Check that `2_process/out/doy_107072210.csv` and `2_process/out/doy_120020150.csv` files were created
- Run `snakemake --cores 1 --delete-all-output` if you want

## Issue 5
- Change Snakefile, combine_site_files.py, and plot_doy_mean.py
- New Snakefile contents:
    * New rule all at the top of the file
```
rule all: 
    input: 
        "3_plot/out/doy_plot.png"
```
    * New rule: combine_site_files
```
rule combine_site_files:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
    output:
        out_file = "2_process/out/combined_doy.csv"
    script:
        "2_process/combine_site_files.py"
```
    * New rule: plot_doy_mean
```
rule plot_doy_mean:
    input:
        in_file = "2_process/out/combined_doy.csv"
    output:
        out_file = "3_plot/out/doy_plot.png"
    params:
        depths = [0, 1, 2, 5]
    script:
        "3_plot/plot_doy_mean.py"
```
- Changes to combine_site_files.py at bottom:
```
if __name__ == '__main__': 
    out_file = snakemake.output['out_file']
    site_files = snakemake.input
    main(out_file, site_files)
```
- Changes to plot_doy_mean.py at bottom:
```
if __name__ == "__main__":
    combined_doy_means = snakemake.input['in_file']
    out_file = snakemake.output['out_file']
    depths = snakemake.params['depths']
    main(in_file, out_file, depths)
```
- Check out issue-5 branch
- Run `snakemake --cores 1 --delete-all-output` 
    * rule all should be present, so no need to specify outputs/rules
- Run `snakemake --cores 1`
- Open `3_plot/out/doy_plot.png` and check that it looks correct
    * One plot pane per site_id
    * One line per depth
- Run `snakemake --cores 1 --delete-all-output` if you want

