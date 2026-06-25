# Smooth HCP resting-state BOLD data with ciftiTools
#
# Applies 5 mm FWHM smoothing (surface and volume) to all subjects that passed
# FD filtering (valid_combined_subjects_FD.csv) across REST1/REST2 x LR/RL.
# Output: /N/project/BayesianBrainMapping/smoothed_bold/sub-{id}/rfMRI_*_smoothed-5mm.dtseries.nii
#
# Run from the project root, or set sourcedir below before sourcing.

if (!exists("dir_project")) {
  # When run standalone: adjust this path to point to the src/ directory
  sourcedir <- "~/Documents/GitHub/BBM-demo/src"
  source(file.path(sourcedir, "setup.R"))
}

# ── Parameters ────────────────────────────────────────────────────────────────
FWHM      <- 5
sessions  <- c("REST1", "REST2")
encodings <- c("LR", "RL")
dir_out   <- file.path(storage_dir, "smoothed_bold")
# Use the midthickness surfaces for smoothing, as they are more biologically relevant.
surfL_path <- file.path(storage_dir, "HCP", "S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii")
surfR_path <- file.path(storage_dir, "HCP", "S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii")

# ── Setup ─────────────────────────────────────────────────────────────────────
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

subject_ids <- read.csv(
  file.path(dir_data, "priors", "filtering", "valid_combined_subjects_FD.csv")
)$subject_id

cat(sprintf("Subjects: %d | Sessions: %s | Encodings: %s | FWHM: %d mm | Threads: %d\n",
            length(subject_ids),
            paste(sessions,  collapse = "/"),
            paste(encodings, collapse = "/"),
            FWHM,
            nThreads))

# ── Smoothing loop (parallelized over subjects) ───────────────────────────────
cl <- makeCluster(nThreads)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

foreach(
  subject = subject_ids,
  .packages = "ciftiTools",
  .export   = c("dir_out", "dir_HCP", "sessions", "encodings", "FWHM", "wb_path",
                "surfL_path", "surfR_path")
) %dopar% {

  ciftiTools.setOption("wb_path", wb_path)

  sub_label   <- paste0("sub-", subject)
  sub_out_dir <- file.path(dir_out, sub_label)
  dir.create(sub_out_dir, recursive = TRUE, showWarnings = FALSE)

  for (session in sessions) {
    for (encoding in encodings) {

      run_tag   <- sprintf("rfMRI_%s_%s", session, encoding)
      bold_path <- file.path(dir_HCP, subject, "MNINonLinear", "Results",
                             run_tag,
                             paste0(run_tag, "_Atlas_MSMAll_hp2000_clean.dtseries.nii"))
      out_path  <- file.path(sub_out_dir,
                             paste0(run_tag, "_Atlas_MSMAll_hp2000_clean",
                                    "_smoothed-", FWHM, "mm.dtseries.nii"))

      if (file.exists(out_path)) next

      if (!file.exists(bold_path)) {
        warning(sprintf("File not found, skipping: %s", bold_path))
        next
      }

      xii <- read_cifti(bold_path)
      xii <- smooth_cifti(xii, surf_FWHM = FWHM, vol_FWHM = FWHM,
                          surfL_fname = surfL_path, surfR_fname = surfR_path)
      write_cifti(xii, out_path)
    }
  }
}

cat("Done.\n")
