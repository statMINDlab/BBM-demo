# Renders the spatial MSE maps from the plot bundle produced by priors_MSE.R.
# Run this on a machine with a working OpenGL device (rgl) -- e.g. a laptop --
# after copying data/outputs/priors_MSE/priors_MSE_bundle.rds over from a
# headless cluster run.

library(ciftiTools)
library(BayesBrainMap)

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

out_dir <- file.path(dir_data, "outputs", "priors_MSE")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- readRDS(file.path(out_dir, "priors_MSE_bundle.rds"))
scaffold <- bundle$scaffold

make_mse_prior <- function(scaffold, sq) {
  mp <- scaffold
  mp$prior <- list(mean = sq)
  class(mp) <- "prior.cifti"
  mp
}

for (combo in bundle$combos) {
  cat(sprintf("Rendering nSubs=%d, smoothing=%s\n", combo$nSubs, combo$smoothing_label))
  mse_prior <- make_mse_prior(scaffold, combo$sq)
  Q <- ncol(mse_prior$prior$mean)
  for (i in seq_len(Q)) {
    label_name <- rownames(mse_prior$template_parc_table)[mse_prior$template_parc_table$Key == i]
    fname <- file.path(plot_dir, paste0(combo$base_name, "_", label_name))
    title <- sprintf("MSE vs original prior | nSubs=%d, smoothing=%s | %s",
                     combo$nSubs, combo$smoothing_label, label_name)
    cat("  Plotting MSE map for component:", label_name, "\n")
    plot(mse_prior, stat = "mean", fname = fname, idx = i, title = title)
  }
}

cat("Done. Spatial maps saved to", plot_dir, "\n")
