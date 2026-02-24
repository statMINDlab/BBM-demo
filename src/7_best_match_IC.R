#Dependencies ##########################################################

#remove.packages("fMRItools")
#devtools::install_github("mandymejia/fMRItools", "7.0")
library(fMRItools)
library("scales")

# Helper functions ################################################

cols <- ciftiTools::ROY_BIG_BL()
jet_diverging <- function(...) {
  scale_fill_gradientn(
    colours = c("cyan", "green", "purple", "blue", "black", "red", "orange", "yellow"),
    values  = rescale(c(-0.5, -0.45, -0.3, -0.15, 0, 0.2, 0.3, 0.5)), 
    limits  = c(-0.5, 0.5),
    oob     = squish
  )
}

# Dice coefficient
dice_coef <- function(a, b) {
  # a and b: binary matrices (V x Q)
  
  if (nrow(a) != nrow(b)) {
    stop("a and b must have the same number of rows (vertices)")
  }
  
  # Ensure numeric 0/1
  a <- as.matrix(a)
  b <- as.matrix(b)
  
  # Intersection counts (Q x Q)
  # crossprod gives t(a) %*% b
  intersect <- crossprod(a, b)
  
  # Column sizes
  size_a <- colSums(a)  # length Q
  size_b <- colSums(b)  # length Q
  
  # Denominator matrix (Q x Q)
  denom <- outer(size_a, size_b, "+")
  
  # Dice matrix
  dice <- ifelse(denom == 0, 0, 2 * intersect / denom)
  
  return(dice)
}

# Threshold maps (binarize)
threshold_maps <- function(mat, thr = 0.5) {
  # threshold absolute value
  bin <- abs(mat) > thr
  storage.mode(bin) <- "numeric"
  bin
}

# Start loading Base Yeo Parcellation #####################################################################

order = list()

# SPATIAL OVERLAP MATRICES - Load base parcellation (YEO 17) ######

# PARCELLATION
parcellation <- readRDS(file.path(dir_data, "outputs", "parcellations", "Yeo17_simplified_mwall.rds"))

v <- c(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)

parcel_ids  <- 1:17

yeo_onehot <- matrix(0, nrow = length(v), ncol = length(parcel_ids))

for (i in seq_along(parcel_ids)) {
  yeo_onehot[, i] <- as.integer(v == parcel_ids[i])
}


# Collapse to 8 networks
mega_map <- list(
  Vis   = c(1,2),
  SomMot   = c(3,4),
  DorsAttn = c(5,6),
  SalVentAttn = c(7,8),
  Limbic   = c(9),
  Cont  = c(10,11,12),
  Default  = c(13,14,15,16),
  TempPar  = c(17)
)

mega_names <- names(mega_map)

yeo_mega <- matrix(0, nrow = nrow(yeo_onehot), ncol = length(mega_map))

for (i in seq_along(mega_map)) {
  yeo_mega[, i] <- rowSums(yeo_onehot[, mega_map[[i]], drop = FALSE]) > 0
}

colnames(yeo_mega) <- mega_names
storage.mode(yeo_mega) <- "numeric"

labs <- rownames(parcellation$meta$cifti$labels$parcels)[parcellation$meta$cifti$labels$parcels$Key > 0]

# Start for MSC #####################################################################

# PARCELLATION
parcellation <- readRDS(file.path(dir_data, "outputs", "parcellations" , "MSC_parcellation.rds"))
basename = "MSC"

v <- c(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)

parcel_ids  <- 1:17

one_hot <- matrix(0, nrow = length(v), ncol = length(parcel_ids))

for (i in seq_along(parcel_ids)) {
  one_hot[, i] <- as.integer(v == parcel_ids[i])
}

dice_mat = dice_coef(one_hot, yeo_mega) # Compute Dice coefficient matrix

best_match <- apply(dice_mat, 1, which.max)

# makes vector with order of indices for plotting
ordered_idx <- order(best_match)

order[[basename]] <- list(ic_order = ordered_idx)


# PARCELLATION
# Start for GICA15 #####################################################################

# PARCELLATION
scalar_input <- read_cifti(file.path(dir_data, "inputs", "GICA15.dscalar.nii"))
basename = "GICA15"

# keep only non-medial wall cortical data
matrix_input <- rbind(scalar_input$data$cortex_left, scalar_input$data$cortex_right)

parcellation <- threshold_maps(matrix_input, thr = colMeans(matrix_input) + 2 * apply(matrix_input, 2, sd))

parcel_ids  <- 1:15

one_hot <- parcellation

dice_mat = dice_coef(one_hot, yeo_mega) # Compute Dice coefficient matrix

best_match <- apply(dice_mat, 1, which.max)

# makes vector with order of indices for plotting
ordered_idx <- order(best_match)

order[[basename]] <- list(ic_order = ordered_idx)

# Start for GICA25 #####################################################################
scalar_input <- read_cifti(file.path(dir_data, "inputs", "GICA25.dscalar.nii"))
basename = "GICA25"

# keep only non-medial wall cortical data
matrix_input <- rbind(scalar_input$data$cortex_left, scalar_input$data$cortex_right)

parcellation <- threshold_maps(matrix_input, thr = colMeans(matrix_input) + 2 * apply(matrix_input, 2, sd))

parcel_ids  <- 1:25

one_hot <- parcellation

dice_mat = dice_coef(one_hot, yeo_mega) # Compute Dice coefficient matrix

best_match <- apply(dice_mat, 1, which.max)

# makes vector with order of indices for plotting
ordered_idx <- order(best_match)

order[[basename]] <- list(ic_order = ordered_idx)

# PROFUMO #####################################################################
profumo <- read_cifti(file.path(dir_data, "inputs", "PROFUMO.dscalar.nii"))
# this prior has the mdial wall not defined for some reason, applying fix.
profumo$meta$cortex$medial_wall_mask <- scalar_input$meta$cortex$medial_wall_mask
basename = "PROFUMO"

# keep only non-medial wall cortical data
matrix_input <- rbind(profumo$data$cortex_left[profumo$meta$cortex$medial_wall_mask$left,], profumo$data$cortex_right[profumo$meta$cortex$medial_wall_mask$right,])

parcellation <- threshold_maps(matrix_input, thr = colMeans(matrix_input) + 2 * apply(matrix_input, 2, sd))

parcel_ids  <- 1:12

one_hot <- parcellation

dice_mat = dice_coef(one_hot, yeo_mega) # Compute Dice coefficient matrix

best_match <- apply(dice_mat, 1, which.max)

# makes vector with order of indices for plotting
ordered_idx <- order(best_match)

order[[basename]] <- list(ic_order = ordered_idx)
