#' shrinkage_summary.R
#' Companion to shrinkage_estimation.r: loads the saved posterior estimates for
#' all subjects, averages the shrinkage (lambda) and subject mean maps over
#' subjects, saves both average maps, and plots shrinkage as a function of the
#' subject mean.
#'
#'

###################
# Set parameters
dir_project = "~/Documents/GitHub/BBM-demo"
sourcedir = file.path(dir_project, "src")
# Run parameter setup
source(file.path(sourcedir, "setup.R"))

smoothing <- 5

# Subjects to sweep over
test_subs <- c("136227", "186444", "275645","729254", "110411", "185341") #"100610", "100408","101006", "101107")

# Sessions to sweep over; both phase-encoding runs (LR and RL) of each
# session are fit jointly, matching 09_fit_BBM.R
sessions <- c("REST1")
encodings <- c("LR", "RL")

# set a base directory to save results and plots; a subdirectory is created
# per variance inflation factor below
base_out_dir <- file.path(dir_project, "manuscript", "shrinkage")

# directory holding the per-subject .rds files written by shrinkage_estimation.r
out_dir <- file.path(base_out_dir, paste0("fromA"))

# where the group averages go
avg_dir <- file.path(out_dir, "average")
dir.create(avg_dir, recursive = TRUE, showWarnings = FALSE)

# vertices are subsampled for the scatterplot only; the averages use all of them
max_points <- 2e5

###################
# Get the prior mean
prior_fname <- if (!is.null(smoothing)) file.path(dir_data, paste0("priors/Yeo17_smooth-", smoothing, "/prior_combined_Yeo17_noGSR_local.rds")) else file.path(dir_data, paste0("priors/Yeo17/prior_combined_Yeo17_noGSR_local.rds"))

prior <- readRDS(prior_fname)

prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

###################
# Average shrinkage and subject mean over subjects
#
# Only the fit_BBM object is saved by shrinkage_estimation.r, so lambda is
# recomputed here exactly as it is defined there:
#   d1 = |MLE - prior mean|, d3 = |MLE - posterior mean|, lambda = d3/d1
# lambda is undefined where the MLE sits on the prior mean (d1 = 0); those
# vertices are dropped per subject rather than poisoning the average, so the
# maps carry a per-vertex count.

lambda_sum <- NULL
lambda_n <- NULL
mean_sum <- NULL
mle_sum <- NULL
n_used <- 0L
scaffold <- NULL

for (test_sub in test_subs) {

  for (session in sessions) {

    ses_label <- session

    fname <- file.path(out_dir, paste0("sub-", test_sub, "_ses-", ses_label, "_smooth-", smoothing, "local_SD.rds"))

    if (!file.exists(fname)) {
      warning("Missing estimate for sub-", test_sub, " ", ses_label, " -- skipping")
      next
    }

    cat("Loading sub-", test_sub, " ", ses_label, "\n", sep = "")
    posterior_estimate <- readRDS(fname)

    # keep the first MLE xifti as the scaffold for writing the average maps out
    if (is.null(scaffold)) scaffold <- posterior_estimate$MLE

    mle <- as.matrix(posterior_estimate$MLE)
    subj_mean <- as.matrix(posterior_estimate$subjNet_mean)

    # distance 1, difference between MLE and prior mean
    d1 = abs(mle - prior$prior$mean)

    # distance 3, difference between MLE and posterior
    d3 = abs(mle - subj_mean)

    # shrinkage is defined as the ratio between distance(MLE, posterior) and distance(prior, MLE)
    lambda = d3/d1
    lambda[!is.finite(lambda)] <- NA

    if (is.null(lambda_sum)) {
      lambda_sum <- matrix(0, nrow(lambda), ncol(lambda))
      lambda_n <- matrix(0L, nrow(lambda), ncol(lambda))
      mean_sum <- matrix(0, nrow(subj_mean), ncol(subj_mean))
      mle_sum <- matrix(0, nrow(mle), ncol(mle))
    }

    lambda_sum <- lambda_sum + ifelse(is.na(lambda), 0, lambda)
    lambda_n <- lambda_n + !is.na(lambda)
    mean_sum <- mean_sum + subj_mean
    mle_sum <- mle_sum + mle

    n_used <- n_used + 1L
  }
}

stopifnot("No saved estimates found in manuscript/shrinkage/fromA" = n_used > 0)

avg_lambda <- lambda_sum / lambda_n   # NaN where no subject had a defined lambda
avg_mean <- mean_sum / n_used
avg_mle <- mle_sum / n_used

cat("Averaged", n_used, "subject-sessions;", nrow(avg_lambda), "greyordinates x",
    ncol(avg_lambda), "components\n")

###################
# Save the average maps

suffix <- paste0("_smooth-", smoothing, "local_SD")

saveRDS(avg_lambda, file.path(avg_dir, paste0("average_shrinkage", suffix, ".rds")))
saveRDS(avg_mean, file.path(avg_dir, paste0("average_mean", suffix, ".rds")))
saveRDS(avg_mle, file.path(avg_dir, paste0("average_MLE", suffix, ".rds")))

avg_lambda_xii <- newdata_xifti(scaffold, ifelse(is.finite(avg_lambda), avg_lambda, NA))
avg_mean_xii <- newdata_xifti(scaffold, avg_mean)

write_cifti(avg_lambda_xii, file.path(avg_dir, paste0("average_shrinkage", suffix, ".dscalar.nii")))
write_cifti(avg_mean_xii, file.path(avg_dir, paste0("average_mean", suffix, ".dscalar.nii")))

# network labels, as used for the component axis of the scatterplot
labels <- vapply(
  seq_len(ncol(avg_lambda)),
  function(i) rownames(prior$template_parc_table)[prior$template_parc_table$Key == i],
  character(1)
)

# same idx / zlim conventions as the per-subject plots in shrinkage_estimation.r
plot(avg_mean_xii, idx = 14,
     zlim = c(-0.3, 0.3),
     title = sprintf("Average posterior mean | %d subjects", n_used),
     fname = file.path(avg_dir, paste0("average_mean", suffix, ".png")))
plot(avg_lambda_xii, idx = 14, zlim = c(0, 1),
     title = sprintf("Average shrinkage (lambda) | %d subjects", n_used),
     fname = file.path(avg_dir, paste0("average_shrinkage", suffix, ".png")))

###################
# Scatterplot: shrinkage as a function of subject mean

scatter <- tibble(
  subject_mean = as.vector(avg_mean),
  shrinkage = as.vector(avg_lambda),
  component = factor(rep(labels, each = nrow(avg_mean)), levels = labels)
)
scatter <- scatter[is.finite(scatter$subject_mean) & is.finite(scatter$shrinkage), ]
saveRDS(scatter, file.path(avg_dir, paste0("shrinkage_vs_mean", suffix, ".rds")))

plot_dat <- scatter
if (nrow(plot_dat) > max_points) plot_dat <- plot_dat[sample(nrow(plot_dat), max_points), ]

p <- ggplot(plot_dat, aes(x = subject_mean, y = shrinkage, color = component)) +
  geom_point(size = 0.3, alpha = 0.2) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    x = "Average subject mean (posterior)",
    y = "Average shrinkage (lambda)",
    color = "Network",
    title = sprintf("Shrinkage vs subject mean | %d subjects, smoothing %dmm", n_used, smoothing)
  ) +
  theme_minimal()

ggsave(file.path(avg_dir, paste0("shrinkage_vs_mean", suffix, ".png")), p,
       width = 8, height = 6, dpi = 300)

########################
# Bins

##Equally-spaced bins for each network
# n_bins <- 30
# 
# trend <- scatter %>%
#   group_by(component) %>%
#   mutate(bin = ntile(subject_mean, n_bins)) %>%
#   group_by(component, bin) %>%
#   summarise(
#     x  = median(subject_mean),
#     lo = quantile(shrinkage, 0.25),
#     hi = quantile(shrinkage, 0.75),
#     y  = median(shrinkage),
#     .groups = "drop"
#   )


## Homogeneous bins across networks
n_bins <- 30

# Common bin edges for ALL components: smallest -> largest of subject_mean
bin_breaks <- seq(
  min(scatter$subject_mean, na.rm = TRUE),
  max(scatter$subject_mean, na.rm = TRUE),
  length.out = n_bins + 1
)

trend <- scatter %>%
  filter(! component %in% c("LimbicA", "LimbicB")) %>%
  mutate(bin = cut(subject_mean, breaks = bin_breaks, include.lowest = TRUE)) %>%
  group_by(component, bin) %>%
  summarise(
    x  = median(subject_mean),
    lo = quantile(shrinkage, 0.25),
    hi = quantile(shrinkage, 0.75),
    y  = median(shrinkage),
    .groups = "drop"
  )

# every network's curve, greyed, repeated behind each panel for comparison
# (no `component` column, so it is not split by the facet)
backdrop <- trend %>% select(x, y, bg = component)

# canonical Yeo17 colours from the parcellation label table when available
network_colors <- if (all(c("Red", "Green", "Blue") %in% names(prior$template_parc_table))) {
  setNames(
    with(prior$template_parc_table, rgb(Red, Green, Blue)),
    rownames(prior$template_parc_table)
  )[labels]
} else {
  setNames(hcl.colors(length(labels), "Dark 3"), labels)
}

p2 <- ggplot(trend, aes(abs(x), y, color = component, fill = component)) +
  geom_line(data = backdrop, aes(x = abs(x), y = y, group = bg), inherit.aes = FALSE,
            color = "grey85", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ component, ncol = 5) +
  scale_color_manual(values = network_colors, guide = "none") +
  scale_fill_manual(values = network_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Average subject mean (posterior)",
    y = "Average shrinkage (lambda)",
    title = sprintf("Shrinkage vs subject mean | %d subjects, smoothing %dmm", n_used, smoothing),
    subtitle = "Median with IQR per bin; grey = all networks"
  ) +
  theme_minimal()

ggsave(file.path(avg_dir, paste0("shrinkage_vs_mean_facet", suffix, ".png")), p2,
       width = 11, height = 7, dpi = 300)

########################
# Same trendlines, all networks overlaid on a single panel
#
# Reuses the bins and colours above. No IQR ribbons here -- 17 overlapping
# bands would be unreadable; the faceted version above carries the spread.

p3 <- ggplot(trend, aes(x, y, color = component)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = network_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Average subject mean (posterior)",
    y = "Average shrinkage (lambda)",
    color = "Network",
    title = sprintf("Shrinkage vs subject mean | %d subjects, smoothing %dmm", n_used, smoothing),
    subtitle = sprintf("Median shrinkage in %d equal-count bins per network", n_bins)
  ) +
  theme_minimal()

ggsave(file.path(avg_dir, paste0("shrinkage_vs_mean_overlay", suffix, ".png")), p3,
       width = 9, height = 6, dpi = 300)

#############################################
#
# Make average plot for reviewers response.
#

# Calculate average within each bin 
average_trend <- trend %>% 
  group_by(bin) %>% 
  summarize(x = mean(x), y = mean(y))

p_condensed <- ggplot() +
  geom_line(data=trend, aes(x, y, group = component, color = component), alpha = 0.5, linewidth = 0.75) +
  geom_line(data = average_trend, aes(x, y), linewidth = 3, color = 'black') +
  scale_color_manual(values = network_colors, guide = "none") +
  scale_fill_manual(values = network_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 0.6), xlim = c(-0.8, 0.8)) +
  labs(
    x = "Average posterior mean",
    y = "Average shrinkage",
    title = NULL,
    subtitle = NULL
  ) +
  theme_minimal(base_size = 16)

ggsave(file.path(avg_dir, paste0("shrinkage_vs_mean_summary", suffix, ".png")), p_condensed,
       width = 11, height = 7, dpi = 300)

cat("Done. Averages and plots saved to", avg_dir, "\n")


