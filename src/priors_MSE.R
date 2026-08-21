# Computes vertex-wise MSE between priors estimated with fewer subjects
# and/or spatial smoothing, relative to the original (348-subject,
# unsmoothed) prior. Outputs a summary table of mean MSE per component
# and spatial MSE maps for each prior x component.
#
# Rendering CIFTI surface maps needs a working OpenGL device (rgl). On
# headless HPC nodes whose rgl was built without OpenGL, this script still
# computes the summary table and saves a lightweight plot bundle; the maps
# can then be rendered on a GL-capable machine via priors_MSE_render.R.

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

# Can this machine render OpenGL surface maps? (FALSE on headless rgl builds)
can_render <- tryCatch(
  requireNamespace("rgl", quietly = TRUE) && !isTRUE(rgl::rgl.useNULL()),
  error = function(e) FALSE
)
if (!can_render) {
  message("No OpenGL device available -- computing table + saving plot bundle only. ",
          "Render maps later with priors_MSE_render.R on a GL-capable machine.")
}

# Minimal prior.cifti carrying the squared-error map in place of the mean,
# reusing the package's plot machinery (plot.prior.cifti) for spatial maps.
make_mse_prior <- function(scaffold, sq) {
  mp <- scaffold
  mp$prior <- list(mean = sq)
  class(mp) <- "prior.cifti"
  mp
}

# Shared plotting scaffold (same template/greyordinates for every prior)
scaffold <- list(
  mask_input = reference$mask_input,
  mask = reference$mask,
  params = reference$params,
  dat_struct = reference$dat_struct,
  template_parc_table = reference$template_parc_table
)

render_maps <- function(mse_prior, base_name, nSubs, smoothing_label) {
  Q <- ncol(mse_prior$prior$mean)
  for (i in seq_len(Q)) {
    label_name <- rownames(mse_prior$template_parc_table)[mse_prior$template_parc_table$Key == i]
    fname <- file.path(plot_dir, paste0(base_name, "_", label_name))
    title <- sprintf("MSE vs original prior | nSubs=%d, smoothing=%s | %s",
                     nSubs, smoothing_label, label_name)
    cat("  Plotting MSE map for component:", label_name, "\n")
    plot(mse_prior, stat = "mean", fname = fname, idx = i, title = title)
  }
}

mse_tbl <- tibble(
  nSubs = integer(),
  smoothing = character(),
  component = character(),
  mean_mse = numeric()
)

# Squared-error maps per prior, kept for the plot bundle
bundle_combos <- list()

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

    base_name <- sprintf("MSE_%s_nSubs-%d_smooth-%s", parcellation, nSubs, smoothing_label)

    # Keep the squared-error map for offline rendering
    bundle_combos[[base_name]] <- list(
      nSubs = nSubs,
      smoothing_label = smoothing_label,
      base_name = base_name,
      sq = sq_err
    )

    if (can_render) {
      render_maps(make_mse_prior(scaffold, sq_err), base_name, nSubs, smoothing_label)
    }
  }
}

# Save summary table
saveRDS(mse_tbl, file.path(out_dir, "priors_MSE_summary.rds"))
write.csv(mse_tbl, file.path(out_dir, "priors_MSE_summary.csv"), row.names = FALSE)

# Save the plot bundle (scaffold + squared-error maps) so the spatial maps can
# be rendered on a GL-capable machine even if this run couldn't render them.
saveRDS(
  list(scaffold = scaffold, combos = bundle_combos),
  file.path(out_dir, "priors_MSE_bundle.rds")
)

cat("Done. Summary table saved to", out_dir, "\n")
if (can_render) {
  cat("Spatial maps saved to", plot_dir, "\n")
} else {
  cat("Plot bundle saved to", file.path(out_dir, "priors_MSE_bundle.rds"),
      "-- render with priors_MSE_render.R\n")
}
