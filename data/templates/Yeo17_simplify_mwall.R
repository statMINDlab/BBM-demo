# Prepare Yeo17 Parcellation for Prior Estimation  

yeo17 <- load_parc("Yeo_17")

# Simplify labels
y <- rownames(yeo17$meta$cifti$labels[[1]])
z <- gsub('17Networks_LH_|17Networks_RH_', '', y)
z <- gsub('_.*', '', z)
z <- factor(z, levels=z[!duplicated(z)])
# Parcellation with new labels
yeo17_simplify <- convert_to_dlabel(
  newdata_xifti(yeo17, as.numeric(z)[c(as.matrix(yeo17)) + 1]),
  levels_old = as.numeric(z)[!duplicated(z)],
  levels = as.numeric(z)[!duplicated(z)] - 1,
  labels = levels(z),
  colors = rgb(
    yeo17$meta$cifti$labels[[1]]$Red,
    yeo17$meta$cifti$labels[[1]]$Green,
    yeo17$meta$cifti$labels[[1]]$Blue,
    yeo17$meta$cifti$labels[[1]]$Alpha
  )[!duplicated(z)],
  add_white = FALSE
)

# Medial wall
yeo17_simplify$meta$cifti$labels[[1]] <- rbind(
  data.frame(Key=-1, Red=1, Green=1, Blue=1, Alpha=0, row.names='BOLD_mwall'),
  yeo17_simplify$meta$cifti$labels[[1]]
)
# Medial wall mask
mwall_path <- file.path(dir_data, "inputs", "Human.MedialWall_Conte69.32k_fs_LR.dlabel.nii")
mwall_cifti <- read_cifti(mwall_path)
mwall_L <- mwall_cifti$data$cortex_left == 0
mwall_R <- mwall_cifti$data$cortex_right == 0
yeo17_simplify$data$cortex_left[!mwall_L,] <- NA
yeo17_simplify$data$cortex_right[!mwall_R,] <- NA

yeo17_simplify_mw <- move_to_mwall(yeo17_simplify, values = c(NA))

saveRDS(yeo17_simplify_mw, file.path(dir_data, "outputs", "Yeo17_simplified_mwall.rds"))

