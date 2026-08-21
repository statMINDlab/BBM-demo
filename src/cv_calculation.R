#' prior_cv_calculation.R
#'

###################
# Set parameters
dir_project = "~/Documents/GitHub/BBM-demo"
sourcedir = file.path(dir_project, "src")
# Run parameter setup
source(file.path(sourcedir, "setup.R"))

smoothing <- 5

# set a base directory to save results and plots; a subdirectory is created
# per variance inflation factor below
base_out_dir <- file.path(dir_project, "manuscript", "cv_calculation")

# directory holding the per-subject .rds files written by shrinkage_estimation.r
out_dir <- file.path(base_out_dir)

###################
# Get the prior mean
prior_fname <- if (!is.null(smoothing)) file.path(dir_data, paste0("priors/Yeo17_smooth-", smoothing, "/prior_combined_Yeo17_noGSR_local.rds")) else file.path(dir_data, paste0("priors/Yeo17/prior_combined_Yeo17_noGSR_local.rds"))

prior <- readRDS(prior_fname)

prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

cv_prior <- sqrt(prior$prior$varNN) / abs(prior$prior$mean)

mean_vec <- as.vector(prior$prior$mean)
sd_vec   <- sqrt(as.vector(prior$prior$varNN))

# build xifti

mean_thr_value <- 0.1

cv_prior[prior$prior$mean < mean_thr_value] <- NA

cv_xifti <- newdata_xifti(prior$dat_struct, cv_prior)

plot(cv_xifti, idx = 14, zlim = c(0.25, 0.5), colors = "plasma" , fname = file.path(out_dir, "cv_plot-plasma.png"))

# thresholded SD
sd_thr <- sqrt(prior$prior$varNN)

sd_thr[prior$prior$mean < mean_thr_value] <- NA

sd_xifti <- newdata_xifti(prior$dat_struct, sd_thr)

plot(sd_xifti, idx = 14, zlim = c(0.06, 0.18), fname = file.path(out_dir, "sd_plot.png"))

# thresholded mean
mean_thr <- prior$prior$mean

mean_thr[prior$prior$mean < mean_thr_value] <- NA

mean_xifti <- newdata_xifti(prior$dat_struct, mean_thr)

plot(mean_xifti, idx = 14, fname = file.path(out_dir, "mean_plot.png"))

##########
# Build data frame

parc_tab <- prior$template_parc_table
parc_tab <- parc_tab[parc_tab$Key > 0, ]
nm <- rownames(parc_tab)

V <- nrow(prior$prior$mean)
stopifnot(length(nm) == ncol(prior$prior$mean))

Q <- ncol(prior$prior$mean)

df <- data.frame(
  mean    = mean_vec,
  sd      = sd_vec,
  cv      = sd_vec / abs(mean_vec),
  network = factor(rep(seq_len(Q), each = V))
)
df$network <- factor(rep(nm, each = V), levels = nm)

df <- df %>%
  filter(mean > 0)

# Colors

network_colors <- setNames(
  with(parc_tab, rgb(Red, Green, Blue)),
  labels
)


#######
# plot
  
p <- ggplot(df, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.1, shape = ".") +
  facet_wrap(~ network, scale = "free") +
  labs(x = "Prior mean", y = "prior SD") + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 10)) +
  theme_minimal()

#p

ggsave(file.path(out_dir, "sd_plot_smooth.png"))
