# Setup

# install.packages("devtools")
# install.packages("gsignal")
# install.packages("ggcorrplot")
# install.packages("ciftiTools")            
# devtools::install_github("mandymejia/fMRIscrub", "14.0")          
# install.packages("fMRItools") # deprecated for new BBM
# devtools::install_github("mandymejia/ciftiTools", "20.0", force=TRUE)
# install.packages("viridis")
# install.packages("BayesBrainMap")
# install.packages("doParallel")
# devtools::install_github("diegoderman/BayesBrainMap", ref = "2.0")

# Load packages
#devtools::load_all("~/Documents/GitHub/BayesBrainMap") # load mean_local normalization version.
#devtools::install("~/Documents/GitHub/BayesBrainMap", quiet = TRUE, upgrade = "never") # install development version of BayesBrainMap
library(fMRItools)       # version 0.8.0
library(ggcorrplot)      # version 0.1.4.1
library(gsignal)         # version 0.3.7
library(ciftiTools)      # version 0.20.0
library(fMRIscrub)       # version 0.14.7
library(viridis)         # version 0.6.5
library(BayesBrainMap)   # version: 0.2.0
library(tidyverse)       # version: 2.0.0
library(purrr)           # version: 0.2.0
library(doParallel)
library(foreach)


# Set CIFTI Workbench path according to the system
if (Sys.info()["sysname"] == "Darwin") {
  wb_path <- "~/workbench/bin_macosxub"
  storage_dir <- "~/Documents/BayesianBrainMapping"
} else if (Sys.info()["sysname"] == "Linux") {
  # if nodename ends with quartz, it's quartz
  if (endsWith(Sys.info()["nodename"], "uits.iu.edu")) {
    wb_path <- "~/Downloads/workbench/bin_rh_linux64"
    storage_dir <- "/N/project/BayesianBrainMapping"
  } else {
    wb_path <- "~/Downloads/workbench/bin_linux64"
    storage_dir <- "~/Documents/BayesianBrainMapping"
  }
} else {
  stop("Unsupported operating system")
}
# Check if the path exists, otherwise throw an error
if (!file.exists(wb_path)) {
  stop(paste("Workbench path does not exist:", wb_path))
}
ciftiTools.setOption("wb_path", wb_path) 

# Set up paths
#########################################################################
#
# INSERT YOUR OWN PATH TO THE HCP DATASET IN THE LINE BELOW TO RUN PIPELINE
#
#########################################################################
# Bring your own HCP access for both restricted and unrestricted data
dir_project <- "~/Documents/GitHub/BBM-demo" # Path to GitHub folder

dir_data <- file.path(dir_project, "data") # Path to data folder

# Bring your own HCP access for both restricted and unrestricted data
# Set CIFTI Workbench path according to the system
if (Sys.info()["sysname"] == "Darwin") {
  dir_HCP <- "~/Documents/hcp_dcwan"
} else if (Sys.info()["sysname"] == "Linux") {
  # if nodename ends with quartz, it's quartz
  if (endsWith(Sys.info()["nodename"], "uits.iu.edu")) {
    dir_HCP <- "/N/project/hcp_dcwan"
    dir_smoothHCP <- "/N/project/BayesianBrainMapping/smoothed_bold"
  } else {
    dir_HCP <- "~/Documents/hcp_dcwan"
  }
} else {
  dir_HCP <- "~/Documents/hcp_dcwan" # set your own directory for HCP
}

HCP_restricted_fname <- file.path(dir_HCP, "..", "restricted_HCP.csv")
# TEST PURPOSES ONLY TRYING WITH RESTRICTED DEMEOGRAPHICS
#HCP_unrestricted_fname <- file.path(dir_HCP, "..", "restricted_HCP_demographics.csv")
HCP_unrestricted_fname <- file.path("~/Documents", "restricted_HCP_demographics.csv")

# Read CSV
#HCP_restricted <- read.csv(HCP_restricted_fname)
HCP_unrestricted <- read.csv(HCP_unrestricted_fname)

# All subject IDS
subject_ids <- HCP_unrestricted$Subject

# Constants
fd_lag_HCP <- 4 # based on multiband factor?
fd_cutoff <- .5 # Motion scrubbing threshold
TR_HCP <- .72 # Repetition time, in seconds
TR_MSC <- 2.2 # Repetition time for MSC data, in seconds
nT_HCP <- 1200 # Timepoints for each resting state scan
min_total_sec <- 600 # Minimum duration of time series after scrubbing (600 sec = 10 min)

# Calculation constants
nThreads = 1 # number of threads to use to estimate priors
Sys.setenv(OMP_NUM_THREADS = "1")

# Parameter sweep definition for prior estimation
encoding_sweep = c("combined") # Using only combined c("LR", "RL", "combined") 
nIC_sweep = c(0, 1, 2, 15, 25) # Yeo, MSC, PROFUMO, GICA 15, GICA 25, see details in 04_estimate_priors.R
GSR_sweep = c(FALSE, TRUE)

# Parameter definition for fit BBM
method_variance = "non-negative"
brainMap_prior = 0 # Yeo 17 selected, see details in 04_estimate_priors.R

