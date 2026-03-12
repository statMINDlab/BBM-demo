# Setup

# install.packages("devtools")
# install.packages("gsignal")
# install.packages("ggcorrplot")
# install.packages("ciftiTools")            
# devtools::install_github("mandymejia/fMRIscrub", "14.0")          
# install.packages("fMRItools") # deprecated for new BBM
# devtools::install_github("mandymejia/fMRItools", "7.0", force=TRUE)
# install.packages("viridis")
# install.packages("BayesBrainMap")
# install.packages("doParallel")
# devtools::install_github("diegoderman/BayesBrainMap", ref = "2.0")

# Load packages
library(ggcorrplot)      # version 0.1.4.1
library(gsignal)         # version 0.3.7
library(ciftiTools)      # version 0.17.4
library(fMRIscrub)       # version 0.14.7
library(viridis)         # version 0.6.5
library(BayesBrainMap)   # version: 0.2.0
library(tidyverse)       # version: 2.0.0
library(purrr)           # version: 0.2.0

# Set CIFTI Workbench path
wb_path <- "/Applications/wb_view.app/Contents/usr/bin"
ciftiTools.setOption("wb_path", wb_path) 

# Set up paths
#########################################################################
#
# INSERT YOUR OWN PATH TO THE HCP DATASET IN THE LINE BELOW TO RUN PIPELINE
#
#########################################################################
# Bring your own HCP access for both restricted and unrestricted data
dir_HCP <- "~/Documents/GitHub/BayesianBrainMapping-priors/data_OSF/inputs/HCP_demo" # Path to folder with HCP demographics CSVs

dir_project <- "~/Documents/GitHub/BayesianBrainMapping-priors" # Path to GitHub folder

dir_data <- file.path(dir_project, "data_OSF") # Path to data folder

# Bring your own HCP access for both restricted and unrestricted data
dir_HCP_demo <- "~/Documents/GitHub/BBM-priors/data_OSF/inputs/HCP_demo" # Path to folder with HCP demographics CSVs

# HCP_unrestricted_fname <- file.path(dir_data, "inputs", "unrestricted_HCP_demographics.csv")
HCP_restricted_fname <- file.path(dir_HCP, "restricted_HCP.csv")
# TEST PURPOSES ONLY TRYING WITH RESTRICTED DEMEOGRAPHICS
HCP_unrestricted_fname <- file.path(dir_data, "inputs", "restricted_HCP_demographics.csv")

# Read CSV
HCP_restricted <- read.csv(HCP_restricted_fname)
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
nThreads = 46 # number of threads to use to estimate priors

# Parameter sweep definition for prior estimation
encoding_sweep = c("combined") # Using only combined c("LR", "RL", "combined") 
nIC_sweep = c(0, 1, 2, 15, 25) # Yeo, MSC, PROFUMO, GICA 15, GICA 25
GSR_sweep = c(FALSE, TRUE)

# Parameter definition for fit BBM
method_variance = "unbiased"
brainMap_prior = 0 # Yeo 17 selected, see details in 5_estimate_priors.R

