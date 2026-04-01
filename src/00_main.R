# Main script to run reproducibility demo in GitHub repository 

# Initialize libraries
Sys.setenv(R_LIBS_USER='$R_LIBS_USER')
###############################################################
#
# CHECK DIRECTORY PATHS IN setup.R BEFORE RUNNING THIS SCRIPT
#
###############################################################
sourcedir = file.path(dir_project, "src")
#sourcedir = "src"

# Setup up dependencies and parameters
source(file.path(sourcedir, "00_setup.R"))

# Run framewise displacement filtering 
source(file.path(sourcedir, "01_fd_time_filtering.R"))

# Filter unrelated subjects
source(file.path(sourcedir, "02_unrelated_filtering.R"))

# Balance sex within age groups
source(file.path(sourcedir, "03_balance_age_sex.R"))

######## Begin estimate priors over the parameter sweep defined in 0_setup.R ######
source(file.path(sourcedir,"04_estimate_prior.R"))

# Intialize performance summary
performance_tbl <- tibble(
  encoding = character(),
  nIC = integer(),
  GSR = logical(),
  elapsed_sec = numeric(),
  user_sec = numeric(),
  sys_sec = numeric(),
  num_thr = numeric()
)

# Go over parameter sweep

for(encoding in encoding_sweep){
  for(nIC in nIC_sweep){
    for(GSR in GSR_sweep){
      
      # run while saving performance data
      timing <- system.time({
        estimate_and_export_prior(encoding,
                                  nIC,
                                  GSR,
                                  dir_data,
                                  TR_HCP,
                                  usePar = nThreads)
      })
      
      performance_tbl <- add_row(
        performance_tbl,
        encoding = encoding,
        nIC = nIC,
        GSR = GSR,
        elapsed_sec = timing["elapsed"],
        user_sec = timing["user.self"],
        sys_sec = timing["sys.self"],
        num_thr = nThreads
      )
    }
  }
}

# save RDS with performance tibble
saveRDS(performance_tbl, file.path(dir_data, "outputs", "prior_estimation_timings.rds"))

##### Begin presenting prior results ######
# visualize prior maps
source(file.path(sourcedir,"05_visualization_prior.R"))

# find best match IC to order FC matrices
source(file.path(sourcedir,"06_best_match_IC.R"))

# visualize FC matrices
source(file.path(sourcedir,"07_visualization_FC.R"))

# visualize matrices of IC overlap
source(file.path(sourcedir,"08_dice_overlap.R"))

##### Begin brain mapping of single HCP subject with Yeo17 priors ########
source(file.path(sourcedir, "09_fit_BBM.R"))
# make visualization of posterior FC estimate maps
source(file.path(sourcedir, "10_BBM_visualization.R"))


