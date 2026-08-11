#' shrinkage_estimation.r
#' Calculate a MLE estimation of the focal subject to estimate the amount of shrinkage
#' caused by the BBM model.
#'
#'

###################
# Set parameters
dir_project = "~/Documents/GitHub/BBM-demo"
sourcedir = file.path(dir_project, "src")
# Run parameter setup
source(file.path(sourcedir, "setup_3.0.R"))

HCP_surf <- read_xifti(surfL_fname = file.path("~/Downloads/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii"),
                                surfR_fname = file.path("~/Downloads/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii"))
smoothing <- 5

# Subjects to sweep over
test_subs <- c("136227", "186444", "275645","729254", "110411", "185341", "100610", "100408","101006", "101107")

# Sessions to sweep over; both phase-encoding runs (LR and RL) of each
# session are fit jointly, matching 09_fit_BBM.R
sessions <- c("REST1")
encodings <- c("LR", "RL")

# set a base directory to save results and plots; a subdirectory is created
# per variance inflation factor below
base_out_dir <- file.path(dir_project, "manuscript", "shrinkage")

###################
# Get the prior mean
prior_fname <- if (!is.null(smoothing)) file.path(dir_data, paste0("priors/Yeo17_smooth-", smoothing, "/prior_combined_Yeo17_noGSR_local.rds")) else file.path(dir_data, paste0("priors/Yeo17/prior_combined_Yeo17_noGSR_local.rds"))

prior <- readRDS(prior_fname)

prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

# Peak concurrency is one worker per (subject, variance inflation) = the two
# nested foreach loops below. On Linux registerDoParallel() uses forked
# workers, which inherit the sourced setup, the loaded prior, and all attached
# packages (same model as 09_fit_BBM.R).
peak_workers <- length(test_subs)
n_cores <- parallel::detectCores()

# distribute threads across each worker's internal fits to avoid
# oversubscription, then register the backend
Sys.setenv(OMP_NUM_THREADS = as.character(4)) #as.character(max(1, floor(n_cores / peak_workers))))
library(doParallel)  # pulls in foreach and parallel

cl <- parallel::makeCluster(peak_workers, type = "FORK")
registerDoParallel(cl)
on.exit(parallel::stopCluster(cl), add = TRUE)

###################
# Simple file logger. Forked workers inherit this function and append progress
# messages (with timestamp and PID) to a shared log file. Watch it live with
# e.g. `tail -f manuscript/shrinkage/run.log`.
dir.create(base_out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(base_out_dir, "run.log")
log_msg <- function(...) {
  cat(sprintf("[%s pid %d] %s\n", format(Sys.time(), "%H:%M:%S"), Sys.getpid(), sprintf(...)),
      file = log_file, append = TRUE)
}

###################
# Outer loop parallelizes over subjects; both LR and RL runs of each session
# are fit jointly (session is an inner loop). The informative posterior does
# not depend on the inflation factor, so it is computed ONCE per subject/session
# here, before the inner variance sweep. Because the inner foreach runs in
# forked workers, they inherit this posterior (and the loaded BOLD) copy-on-write,
# so it is reused across inflation values without recomputation.

lambda_medians <- foreach(
  test_sub = test_subs,
  .combine = "c",
  .packages = c("ciftiTools", "BayesBrainMap"),
  .errorhandling = "pass"
) %dopar% {
  
  ciftiTools.setOption("wb_path", wb_path)   # workers don't inherit this
  log_msg("START subject %s", test_sub)
  
  sapply(sessions, function(session) {
    ses_label <- session
    
    BOLD_fnames <- file.path(
      dir_smoothHCP,
      paste0("sub-", test_sub),
      sprintf("rfMRI_%s_%s_Atlas_MSMAll_hp2000_clean_smoothed-%dmm.dtseries.nii",
              session, encodings, smoothing)
    )
    BOLD_focal <- lapply(BOLD_fnames, read_cifti)
    
    posterior_estimate <- fit_BBM(
      BOLD_focal, prior = prior, var_method = method_variance,
      TR = TR_HCP, drop_first = 5, GSR = FALSE,
      scale_sm_FWHM = 0, usePar = FALSE, MLE = TRUE
    )
    
    out_dir <- file.path(base_out_dir, "fromA")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(posterior_estimate, file.path(out_dir, paste0(
      "sub-", test_sub, "_ses-", ses_label, "_smooth-", smoothing, "local_SD.rds")))
    
    d1 <- abs(as.matrix(posterior_estimate$MLE) - prior$prior$mean)
    d3 <- abs(as.matrix(posterior_estimate$MLE) -
                as.matrix(posterior_estimate$subjNet_mean))
    lambda <- d3 / d1
    
    median(lambda[, 14])
  })
}
