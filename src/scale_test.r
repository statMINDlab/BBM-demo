# Testing new smoothing 

library("BayesBrainMap")
find.package("BayesBrainMap")

dir_HCP <- "/N/project/hcp_dcwan"
dir_smoothHCP <- "h" 
HCP_surf <- read_xifti(surfL_fname = file.path("~/Downloads/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii"),
                       surfR_fname = file.path("~/Downloads/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii"))
smoothing <- 5
test_sub <- "133827"
BOLD_fname <- file.path(dir_HCP, 
                        test_sub, 
                        sprintf("MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii"))

scale_BOLD_fname <- file.path(dir_smoothHCP, 
                         paste0("sub-", test_sub), 
                         sprintf("rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean_smoothed-%dmm.dtseries.nii", smoothing))


# GICA
dir_data <- "~/Documents/GitHub/BBM-demo/data"
GICA <- readRDS(file.path(dir_data, "templates", "Yeo17_simplified_mwall.rds"))

# read BOLD
BOLD <- read_cifti(BOLD_fname, brainstructures=c("left", "right"))

#normalize BOLD

normalized_BOLD <- norm_BOLD(as.matrix(BOLD), TR=TR_HCP, scale="local", hpf=0.01)

normalized_BOLD_cifti <- newdata_xifti(BOLD, normalized_BOLD)

# Apply dual regression

DR_xifti <- dual_reg2(BOLD, template=GICA, scale = 'local')

# plot  sample timeseries

###########################################
#
# New test, prior results
#
###########################################

# Testing smoothing effects on priors.

prior_smooth <- "/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local.rds"
prior_nosmooth <- "/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25/prior_combined_Yeo17_noGSR_local.rds"

smooth <- readRDS(prior_smooth)
nosmooth <- readRDS(prior_nosmooth)

smooth_xii <- newdata_xifti(smooth$dat_struct, smooth$prior$mean)
nosmooth_xii <- newdata_xifti(nosmooth$dat_struct, nosmooth$prior$mean)

for (i in 1:17){
  plot(smooth_xii, zlim = c(-0.5, 0.5), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-5_", i))
  plot(nosmooth_xii, zlim = c(-0.5, 0.5), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-NULL_", i))
  
}

smooth_xii <- newdata_xifti(smooth$dat_struct, sqrt(smooth$prior$varNN))
nosmooth_xii <- newdata_xifti(nosmooth$dat_struct, sqrt(nosmooth$prior$varNN))

for (i in 1:17){
  plot(smooth_xii, zlim = c(0, 0.25), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-5_sd", i))
  plot(nosmooth_xii, zlim = c(0, 0.25), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-NULL_sd", i))
  
}

cv_smooth <- sqrt(smooth$prior$varNN) / as.matrix(smooth$prior$mean)

cv_smooth_xii <- newdata_xifti(smooth$dat_struct, cv_smooth)

cv_nosmooth <- sqrt(as.matrix(nosmooth$prior$varNN)) / as.matrix(nosmooth$prior$mean)

cv_nosmooth_xii <- newdata_xifti(nosmooth$dat_struct, cv_nosmooth)

i = 14
plot(cv_smooth_xii, zlim = c(-10, 10), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-5_cv", i))
plot(cv_nosmooth_xii, zlim = c(-10, 10), idx = i, fname = paste0("/N/u/dderman/Quartz/Documents/GitHub/BBM-demo/data/priors_mean-scale/Yeo17/nSubs-25_smooth-5/prior_combined_Yeo17_noGSR_local_smooth-NULL_cv", i))

summary(abs(cv_local_xii$data$cortex_left[,14]) - abs(cv_robust_xii$data$cortex_left[,14]))

plot(cv_local_xii - cv_robust_xii, idx = 14)