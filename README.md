# BBM-priors

This repository holds population-derived priors/templates for Bayesian Brain Mapping (BBM), also known as template ICA. 

It accompanies the manuscript **Bayesian brain mapping: population-informed individualized functional topography and connectivity**. Nohelia Da Silva Sanchez, Diego Derman, Saige Rutherford, Damon D. Pham, Ellyn R. Butler, Mary Beth Nebel, Amanda F. Mejia
[arXiv preprint](https://doi.org/10.48550/arXiv.2602.01551)

Population-derived priors and supporting data used in this study are publicly available through the Open Science Framework (OSF): https://osf.io/k6vx8/files/osfstorage

## Vignette/Demo<a name="demo-header"></a>

- Step 1: [**Population priors**](https://htmlpreview.github.io/?https://github.com/mandymejia/BBM-priors/refs/heads/main/demo/demo_step1.html). Usage of `estimate_prior()` and manipulations of priors.

- Step 2: [**Subject-level BBM fit**](https://htmlpreview.github.io/?https://github.com/mandymejia/BBM-priors/refs/heads/main/demo/demo_step2.html). Usage of `fit_BBM()`, `id_engagements()`, and manipulation of subject-level maps and results.


## Key Dependencies

- **R packages**: `BayesBrainMap`, `ciftiTools`, `fMRIscrub`, `kableExtra`, `magrittr`, `knitr`, `rgl`, `bookdown`
- **External tools**: Connectome Workbench (`wb_command`) must be installed and path configured via `wb_path`
- **Data format**: CIFTI neuroimaging format (`.dtseries.nii`, `.dlabel.nii`)

## Repository Structure

- `src/` - Contains all the scripts necessary to reproduce and customize the pipeline.
- `src/manuscript/` - Figure generation scripts for the paper.

- `demo/` - Tutorial/demos showcasing different aspects of the pipeline. They can be modified to customize and reproduce the pipeline step-by-step.

- `data/` - Directory for input/output data. When OSF data is required, it will be downloaded here.

## Running the Analysis

To make the package more approachable, we provide different ways of approaching this repository, based on the level of depth and customabilization you need.

**1. View Demo** The [generated htmls](#demo-header) provide a complete overview of each stage of the pipeline

**2. Reproduce Demo analsysis:**

The demo can be fully reproduced without having to estimate population priors and changed. To accomplish this, Yeo17 priors and subject-level results for a single HCP subject are programmaticaly downloaded from OSF using the `osfr` R package.

```r
# Install osfr library
install.packages("osfr")
```

Then, individual snippets can be evaluated. This can be achieved by toggling each snippet to `eval=TRUE`, graphically on RStudio.

**3. Full pipeline**

To run the full pipeline, including prior estimation, you will need access to the HCP and the templates available in `data/templates`. Alternative, you can download and see the pre-estimated priors from the Open Science Foundation [repository](https://osf.io/k6vx8/files/osfstorage). 

```r
# BEFORE RUNNING: you will need to bring your own HCP access, modify dir_hcp path in src/setup.R accordingly.
source("src/00_main.R")
```

Without HCP access, you can still obtain the estimated priors and reproduce most article results. Estimated templates are publicly available in `data_OSF`. To run the pipeline, the relevant `data_OSF` material will be downloaded to the `data` folder. All raw outputs will be written within the `data` folder.

# General BayesBrainMap Package Usage

The demo uses two main functions from `BayesBrainMap`:

1. **`estimate_prior()`** - Estimates group-level statistical priors from training data
   - Requires BOLD fMRI paths, template parcellation, and preprocessing parameters
   - Resource-intensive (~27 hours, 135GB memory for ~350 subjects)

2. **`fit_BBM()`** - Fits subject-level model using pre-estimated priors
   - Input: BOLD data paths and prior object
   - Output: Individualized functional network maps

3. **`id_engagements()`** - Identifies regions of significant deviation from prior mean

## Data Requirements

- HCP resting-state fMRI data in CIFTI format
- HCP demographic data (unrestricted and restricted CSVs) for subject filtering
- Yeo17 or GICA parcellations as template files

## Code Conventions

- Path variables are typically set at script top: `dir_project`, `dir_data`, `dir_HCP`, `wb_path`
- Subject filtering pipeline: FD motion filtering → unrelated filtering → sex/age balancing
- Prior naming convention: `prior_<encoding>_<parcellation>_<GSR>.rds` (e.g., `prior_combined_Yeo17_noGSR.rds`)
