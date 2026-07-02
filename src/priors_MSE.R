# Computes vertex-wise MSE between priors estimated with fewer subjects
# and/or spatial smoothing, relative to the original (348-subject,
# unsmoothed) prior. Outputs a summary table of mean MSE per component
# and spatial MSE maps for each prior x component.

# This script only reads previously-saved prior .rds files, so it loads its
# own minimal dependencies rather than sourcing the estimation setup.R (which
# requires restricted HCP subject data not needed here).
library(ciftiTools)
library(BayesBrainMap)
library(tibble)

dir_project <- "~/Documents/GitHub/BBM-demo"
dir_data <- file.path(dir_project, "data")

# CIFTI Workbench path (mirrors the detection logic in setup.R)
if (Sys.info()["sysname"] == "Darwin") {
  wb_path <- "~/workbench/bin_macosxub"
} else if (endsWith(Sys.info()["nodename"], "uits.iu.edu")) {
  wb_path <- "~/Downloads/workbench/bin_rh_linux64"
} else {
  wb_path <- "~/Downloads/workbench/bin_linux64"
}
ciftiTools.setOption("wb_path", wb_path)

parcellation <- "Yeo17"
encoding <- "combined"
gsr_status <- "noGSR"

nSubs_sweep <- c(25, 50, 100, 200)
smoothing_sweep <- list(NULL, 5)

prior_filename <- sprintf("prior_%s_%s_%s.rds", encoding, parcellation, gsr_status)

# Load reference (original, all 348 subjects, no smoothing) prior
reference_path <- file.path(dir_data, "priors", parcellation, prior_filename)
reference <- readRDS(reference_path)
reference$template_parc_table <- subset(reference$template_parc_table, reference$template_parc_table$Key > 0)
reference_mean <- reference$prior$mean

# Output locations
out_dir <- file.path(dir_data, "outputs", "priors_MSE")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

mse_tbl <- tibble(
  nSubs = integer(),
  smoothing = character(),
  component = character(),
  mean_mse = numeric()
)

for (nSubs in nSubs_sweep) {
  for (smoothing in smoothing_sweep) {

    dir_name <- sprintf("nSubs-%d", nSubs)
    if (!is.null(smoothing)) dir_name <- paste0(dir_name, "_smooth-", smoothing)
    smoothing_label <- if (is.null(smoothing)) "none" else as.character(smoothing)

    prior_path <- file.path(dir_data, "priors", parcellation, dir_name, prior_filename)
    if (!file.exists(prior_path)) {
      warning(sprintf("Missing prior file for nSubs=%d, smoothing=%s (%s) -- skipping",
                       nSubs, smoothing_label, prior_path))
      next
    }

    cat(sprintf("Processing nSubs=%d, smoothing=%s\n", nSubs, smoothing_label))

    prior <- readRDS(prior_path)
    prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

    if (!identical(dim(prior$prior$mean), dim(reference_mean))) {
      warning(sprintf("Dimension mismatch for nSubs=%d, smoothing=%s -- skipping",
                       nSubs, smoothing_label))
      next
    }

    # Vertex-wise squared error (V x Q) against the original prior
    sq_err <- (prior$prior$mean - reference_mean)^2
    component_mse <- colMeans(sq_err)
    Q <- length(component_mse)

    for (i in seq_len(Q)) {
      label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]
      mse_tbl <- add_row(
        mse_tbl,
        nSubs = nSubs,
        smoothing = smoothing_label,
        component = label_name,
        mean_mse = component_mse[i]
      )
    }

    # Reuse the package's CIFTI plotting machinery by substituting the
    # squared-error matrix in place of the prior mean
    mse_prior <- prior
    mse_prior$prior$mean <- sq_err

    base_name <- sprintf("MSE_%s_nSubs-%d_smooth-%s", parcellation, nSubs, smoothing_label)

    for (i in seq_len(Q)) {
      label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]
      fname <- file.path(plot_dir, paste0(base_name, "_", label_name))
      title <- sprintf("MSE vs original prior | nSubs=%d, smoothing=%s | %s",
                        nSubs, smoothing_label, label_name)

      cat("  Plotting MSE map for component:", label_name, "\n")

      plot(
        mse_prior,
        stat = "mean",
        fname = fname,
        idx = i,
        title = title
      )
    }
  }
}

# Save summary table
saveRDS(mse_tbl, file.path(out_dir, "priors_MSE_summary.rds"))
write.csv(mse_tbl, file.path(out_dir, "priors_MSE_summary.csv"), row.names = FALSE)

cat("Done. Summary table and plots saved to", out_dir, "\n")
