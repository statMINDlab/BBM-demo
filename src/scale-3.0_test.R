#' This is a test script to solve normalization problems within BBM.

### Start setting up of environment ############################################

# set source directory
sourcedir = "~/Documents/GitHub/BBM-demo/src"
# run setup steps
# Setup up dependencies and parameters
source(file.path(sourcedir, "setup_3.0.R"))

# source brainMap for function
#source(file.path(sourcedir, "11_brainMap.R"))

### Start parameter definition #################################################

# set output directory
manuscript_output_dir <- "~/Documents/GitHub/BBM-demo/manuscript"
output_dir <- file.path(manuscript_output_dir, "outputs", "brain_map")
# create output directory if it does not exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# set subject and session
subject_ids <- c("100408", "186444") # example subject
session_id <- c("REST1")

# set number of concurrent openMP cores
Sys.setenv(OMP_NUM_THREADS = as.character(48/length(subject_ids)))

# set parameters
encoding <- c("LR", "RL") 
smoothing <- 5 # in mm FWHM
scrubbing <- TRUE
# Define prior path based on selected nIC
nIC <- brainMap_prior
prior_path <- if (nIC == 0) {
  file.path(dir_data, "priors", "Yeo17_smooth-5", "prior_combined_Yeo17_noGSR_local.rds")
} else if (nIC == 1) {
  file.path(dir_data, "priors", "MSC/nSubs-25_smooth-5", "prior_combined_MSC_noGSR.rds")
} else if (nIC == 2) {
  file.path(dir_data, "priors", "PROFUMO", "prior_combined_PROFUMO_noGSR.rds")
} else {
  file.path(dir_data, "priors", sprintf("GICA%d", nIC), paste0("prior_combined_", sprintf("GICA%d", nIC), "_noGSR.rds"))
}

### 1. Check BOLD as it goes into estimate_prior() ##############

BOLD_paths1 <- file.path(dir_smoothHCP, 
                         paste0("sub-", subject_ids), 
                         sprintf("rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean_smoothed-%dmm.dtseries.nii", smoothing))
encoding1 = "LR"
session1 = "REST1"

BOLD_paths2 <- file.path(dir_smoothHCP, 
                         paste0("sub-", subject_ids), 
                         sprintf("rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean_smoothed-%dmm.dtseries.nii", smoothing))
encoding2 = "RL"
session2 = "REST1"

# Start scrubbing procedure, keeping ten minutes below FD threshold.

# get fd flags from inputs
fd_flags = readRDS(file.path(dir_data, "priors", "filtering", "fd_flags.rds"))

# Define number of volumes to keep
keep_volumes <-floor(min_total_sec / TR_HCP)

# make nested list fd_flags into tibble, for easy vectorization
fd_tbl <- fd_flags %>%
  imap_dfr(function(enc_list, encoding) {
    enc_list %>%
      imap_dfr(function(subj_list, subject_id) {
        subj_list %>%
          imap_dfr(function(fd_vec, session) {
            tibble(
              encoding = encoding,
              subject  = subject_id,
              session  = session,
              fd       = list(fd_vec)
            )
          })
      })
  })

# filtered tibble for final subject list 
fd_tbl <- fd_tbl %>%
  filter(subject %in% subject_ids)

# Obtain the indices to keep, all subjects represented by 10 minutes
fd_tbl <- fd_tbl %>%
  mutate(
    fd_scrubbed = map(
      fd,
      ~ {
        over_threshold <- !.x
        (cumsum(over_threshold) <= keep_volumes) & over_threshold
      }
    )
  )

# make sure that I am keeping the required amount of volumes for each subject
stopifnot(all(unlist(lapply(fd_tbl$fd_scrubbed, sum)) == keep_volumes))

# format scrub indices in BBM-friendly way
scrub_BOLD1 <- fd_tbl %>%
  filter(encoding == encoding1 & session == session1) %>%
  select(fd_scrubbed)
scrub_BOLD2 <- fd_tbl %>%
  filter(encoding == encoding2 & session == session2) %>%
  select(fd_scrubbed)
scrub <- list(unlist(scrub_BOLD1, recursive = FALSE), 
              unlist(scrub_BOLD2, recursive = FALSE))


# FD scrubbing

# Yeo17 parcellation

GICA <- readRDS(file.path(dir_data, "templates", "Yeo17_simplified_mwall.rds"))

# Include certain ICs (1:17 not 0 or -1 -> medial wall)
valid_keys <- GICA$meta$cifti$labels[[1]]$Key
inds <- valid_keys[valid_keys > 0]


# Diego Debug print libpaths
scale_parameter = "local"
smoothing = 5
scale_sm_FWHM = 0
GSR = FALSE
usePar = FALSE

prior <- estimate_prior(
  BOLD = BOLD_paths1,
  BOLD2 = BOLD_paths2,
  FC = FALSE,
  scale = scale_parameter, # Added mean local scaling of the input BOLD.
  scale_sm_FWHM = scale_sm_FWHM, # Preprocessed data is assumed smoothed.
  template = GICA,
  GSR = GSR,
  TR = TR_HCP,
  hpf = 0.01,
  Q2 = 0,
  Q2_max = NULL,
  verbose = TRUE,
  inds = inds,
  brainstructures = c("left", "right"),
  drop_first = 15,
  scrub = NULL, #scrub,
  usePar=usePar,
  wb_path = wb_path
)

summary(prior$prior$mean)
plot(newdata_xifti(GICA, prior$prior$mean), idx = 14)

#### Check estimated priors ############################

mean_scaling_fname <- "~/Documents/GitHub/BBM-demo/data/priors/Yeo17_smooth-5/prior_combined_Yeo17_noGSR_local.rds"
sd_scaling_fname <- "~/Documents/GitHub/BBM-demo/data/priors/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR.rds"
gsr_fname <- "~/Documents/GitHub/BBM-demo/data/priors/Yeo17_nSubs-5/prior_combined_Yeo17_GSR.rds"

# plot mean scaled prior
mean_scale_prior <- readRDS(mean_scaling_fname)
sd_scale_prior <- readRDS(sd_scaling_fname)
gsr_prior <- readRDS(gsr_fname)

ms_m <- .2
ms_sd <- .3

plot(newdata_xifti(GICA, mean_scale_prior$prior$mean), zlim = c(-ms_m, ms_m), idx = 14, fname = "prior-Yeo17_nSubs-25_Smooth-5_local.png")
plot(newdata_xifti(GICA, gsr_prior$prior$mean), zlim = c(-ms_m, ms_m), idx = 14, fname = "prior-Yeo17_nSubs-25_Smooth-5_local_GSR.png")

plot(newdata_xifti(GICA, sqrt(mean_scale_prior$prior$varNN)), zlim = c(0, ms_m/2), idx = 14, fname = "prior-Yeo17_nSubs-25_Smooth-5_mean-SD.png")
plot(newdata_xifti(GICA, sqrt(sd_scale_prior$prior$varNN)), zlim = c(0, ms_sd/2), idx = 14, fname = "prior-Yeo17_nSubs-25_Smooth-5_SD-NN.png")

### Coefficient of Variation

plot(newdata_xifti(GICA, sqrt(sd_scale_prior$prior$varNN)/abs(sd_scale_prior$prior$mean)), zlim = c(0, 1), idx = 14, fname = "scale_sd_CV.png")
plot(newdata_xifti(GICA, sqrt(mean_scale_prior$prior$varNN)/abs(mean_scale_prior$prior$mean)), zlim = c(0, 1), idx = 14, fname = "scale_mean_CV.png")

cv_difference <- sqrt(sd_scale_prior$prior$varNN)/abs(sd_scale_prior$prior$mean) - sqrt(mean_scale_prior$prior$varNN)/abs(mean_scale_prior$prior$mean)

plot(newdata_xifti(GICA, cv_difference), zlim = c(-1, 1), idx = 14, fname = "difference_CV.png")


BOLD <- as.matrix(read_xifti(BOLD_paths1))
summary(rowMeans(BOLD))
summary(sqrt(rowVars(BOLD)))
