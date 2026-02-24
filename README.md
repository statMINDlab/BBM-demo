# BayesianBrainMapping-priors

This repository holds population-derived priors/templates for Bayesian Brain Mapping (BBM), also known as template ICA. 

It accompanies the manuscript **Bayesian brain mapping: population-informed individualized functional topography and connectivity**. Nohelia Da Silva Sanchez, Diego Derman, Saige Rutherford, Damon D. Pham, Ellyn R. Butler, Mary Beth Nebel, Amanda F. Mejia
[arXiv preprint](https://doi.org/10.48550/arXiv.2602.01551)

## Key Dependencies

- **R packages**: `BayesBrainMap`, `ciftiTools`, `fMRIscrub`, `kableExtra`, `magrittr`, `knitr`, `rgl`, `bookdown`
- **External tools**: Connectome Workbench (`wb_command`) must be installed and path configured via `wb_path`
- **Data format**: CIFTI neuroimaging format (`.dtseries.nii`, `.dlabel.nii`)

## Repository Structure

- `src/00_main.R` - Main entry point that runs the full reproducibility pipeline
- `src/0_setup.R` - Configuration: package loading, paths, and analysis parameters
- `src/1_fd_time_filtering.R` through `src/6_visualization_prior.R` - Pipeline scripts (run in order)
- `src/manuscript/` - Figure generation scripts for the paper
- `demo/BBM-demo.Rmd` - Tutorial/demo showcasing the BBM workflow
- `data/` - Directory for input/output data 

## Running the Analysis

**Full pipeline:**

To run the full pipeline, including prior estimation, you will need access to the HCP.

```r
# BEFORE RUNNING: you will need to bring your own HCP access, modify dir_hcp path in src/0_setup.R accordingly.
source("src/00_main.R")
```

To reproduce all results after prior estimations, all necessary inputs are publicly available in the Open Science Foundation [repository](https://osf.io/k6vx8). To run the pipeline,`data_OSF` is expected in the root directory. All raw outputs will be written within the `outputs` subfolder.

**Demo only:**

```bash
git clone git@github.com:mandymejia/BBM-priors.git 
```

```r
# Render the demo
rmarkdown::render("demo/BBM-demo.Rmd")
```

## BayesBrainMap Package Usage

The demo uses two main functions from `BayesBrainMap`:

1. **`estimate_prior()`** - Estimates group-level statistical priors from training data
   - Requires BOLD fMRI paths, template parcellation, and preprocessing parameters
   - Resource-intensive (~27 hours, 135GB memory for ~350 subjects)

2. **`BrainMap()`** - Fits subject-level model using pre-estimated priors
   - Input: BOLD data paths and prior object
   - Output: Individualized functional network maps

3. **`engagements()`** - Identifies regions of significant deviation from prior mean

## Data Requirements

- HCP resting-state fMRI data in CIFTI format
- HCP demographic data (unrestricted and restricted CSVs) for subject filtering
- Yeo17 or GICA parcellations as template files

## Code Conventions

- Path variables are typically set at script top: `dir_project`, `dir_data`, `dir_HCP`, `wb_path`
- Subject filtering pipeline: FD motion filtering → unrelated filtering → sex/age balancing
- Prior naming convention: `prior_<encoding>_<parcellation>_<GSR>.rds` (e.g., `prior_combined_Yeo17_noGSR.rds`)
