# Bayesian Brain Mapping demo & priors

Demo for Bayesian Brain Mapping (BBM) to accompany the manuscript:

**Bayesian brain mapping: population-informed individualized functional topography and connectivity**  
Nohelia Da Silva Sanchez, Diego Derman, Saige Rutherford, Damon D. Pham, Ellyn R. Butler, Mary Beth Nebel, Amanda F. Mejia  
Preprint: https://doi.org/10.48550/arXiv.2602.01551

Public data assets (priors, example posteriors):  
https://osf.io/k6vx8/files/osfstorage

## Start Here: How to use this repository

### Option 1: I only want to read examples and understand the workflow
Open the rendered demos (no local run required):

- Demo 1: Population priors  
   https://htmlpreview.github.io/?https://github.com/mandymejia/BBM-priors/refs/heads/main/demo/demo_step1.html
- Demo 2: Subject-level BBM fit  
   https://htmlpreview.github.io/?https://github.com/mandymejia/BBM-priors/refs/heads/main/demo/demo_step2.html

This path is best if you want conceptual understanding of `estimate_prior()`, `fit_BBM()`, and `id_engagements()`.

### Option 2: I want to run the demo step by step using different templates and making custom plots
Use `demo/demo_step1.Rmd` and/or `demo/demo_step2.Rmd` interactively:

1. Install core dependencies (see [Requirements](#requirements)).
2. Install `osfr` to download precomputed files (priors and model fits) from OSF.
3. (Optional) Change template settings (for example `template <- "Yeo17"` to other available templates).
4. Run chunks sequentially.

Important notes:
- Demo code is meant to be educational by default; many computationally heavy chunks or those that depend on HCP data are marked `eval=FALSE`. By default, knitting the demos will use the Yeo17 template and the plots saved within this repo.
- You can run most of the demo without estimating priors from scratch by downloading priors and posteriors from OSF.

### Option 3: I want to run and modify the full pipeline on my own dataset
Use the `src/` scripts. This path is for users doing major customization (new templates, different subject selection, alternative preprocessing, custom outputs).

Recommended order:

1. Configure environment and paths in `src/setup.R`.
2. Run the full chain with `src/00_main.R` or run scripts individually:
 - `01_fd_time_filtering.R`: motion and scan-duration filtering
 - `02_unrelated_filtering.R`: keep unrelated participants
 - `03_balance_age_sex.R`: age/sex balancing
 - `04_estimate_prior.R`: estimate priors over parameter sweep
 - `05_visualization_prior.R` to `08_dice_overlap.R`: prior diagnostics and visualizations
 - `09_fit_BBM.R` and `10_BBM_visualization.R`: subject-level mapping and plots

## Requirements

### Software
- R (tested with modern R 4.x setups)
- Connectome Workbench (`wb_command`) installed and accessible

### R packages
At minimum, the pipeline uses:
- `BayesBrainMap`
- `ciftiTools`
- `fMRIscrub`
- `tidyverse`
- `ggcorrplot`
- `gsignal`
- `viridis`
- `kableExtra`, `magrittr`, `knitr`, `bookdown` (for demos)

For workflows dependent on precalculated priors or posteriors:
- `osfr`

### Data format
- CIFTI inputs (`.dtseries.nii`, `.dlabel.nii`) for HCP-style workflows
- RDS inputs for preformated templates
- Demographics CSV files for participant filtering

## General BayesBrainMap Package Usage

The demo uses two main functions from BayesBrainMap:

1. `estimate_prior()` - Estimates group-level statistical priors from training data
        Requires BOLD fMRI paths, template parcellation, and preprocessing parameters
        Resource-intensive (~27 hours, 135GB memory for ~350 subjects)
2. `fit_BBM()` - Fits subject-level model using pre-estimated priors
        Input: BOLD data paths and prior object
        Output: Individualized functional network maps
3. `id_engagements()` - Identifies regions of significant deviation from prior mean


## Repository Layout

```text
BBM-priors/
├── demo/                 # Tutorial notebooks and rendered HTMLs
├── src/                  # Full reproducible analysis pipeline
│   └── manuscript/       # Figure scripts
├── data/
│   ├── templates/        # Template parcellations/maps
│   ├── priors/           # Estimated or downloaded priors
│   └── posteriors/       # Subject-level BBM outputs
└── manuscript/           # Manuscript figures and outputs
```

## Typical Inputs and Outputs

### Inputs you must define
- Project/data paths (`dir_project`, `dir_data`, `dir_HCP`) in `src/setup.R`
- HCP restricted/unrestricted demographic CSV locations
- Workbench path (`wb_path`)
- Parameter sweeps (`encoding_sweep`, `nIC_sweep`, `GSR_sweep`)

### Main outputs
All outputs are set to `dir_data`:
- Filtering tables in `data/priors/filtering/`
- Priors in `data/priors/<template>/prior_<encoding>_<template>_<GSR|noGSR>.rds`
- Prior plots in `data/priors/<template>/plots_maps/` and `data/priors/<template>/plots_FC/`
- Subject-level outputs in `data/posteriors/sub-<id>/<template>/`

## Customization Tips

- To compare template families (Yeo17, GICA, MSC, PROFUMO), update template/IC settings in `src/setup.R` and rerun from prior estimation onward.
- For faster development iterations, reduce subject lists and disable large parameter sweeps first.
- Keep raw data external; write outputs under `data/` to preserve script compatibility.

## Citation

If you use this repository, please cite:

Da Silva Sanchez N, Derman D, Rutherford S, Pham DD, Butler ER, Nebel MB, Mejia AF.  
*Bayesian brain mapping: population-informed individualized functional topography and connectivity.*  
arXiv (2026). https://doi.org/10.48550/arXiv.2602.01551
