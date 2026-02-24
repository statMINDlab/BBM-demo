dice_overlap <- function(X, Y = NULL) {
    # If Y is NULL: X must be a binary matrix (rows=voxels, cols=networks)
    # If Y is provided: X and Y must be binary vectors (2-network case)

    if (is.null(Y)) {
        if (!is.matrix(X)) stop("X must be a matrix")
        ux <- unique(as.vector(X))
        if (!all(ux %in% c(0, 1))) stop("Matrix X must be binary")

        XtX   <- crossprod(X)             
        sizes <- diag(XtX)
        denom <- outer(sizes, sizes, "+")
        D <- (2 * XtX) / denom

        return(D)

    } else {
        x <- as.vector(X)
        y <- as.vector(Y)
        if (length(x) != length(y)) stop("x and y must have the same length")
        ux <- unique(c(x, y))
        if (!all(ux %in% c(0, 1))) stop("x and y must be binary")

        inter <- sum(x * y)
        ax <- sum(x)                
        by <- sum(y)                
        D <- (2 * inter) / (ax + by)

        return(D)
    }
}

source(file.path(dir_project, "src", "12_best_match_IC.R"))

########## PARCELLATIONS ##########

lim <- 1

# Yeo17
parcellation <- readRDS(file.path(dir_data, "outputs", "parcellations", "Yeo17_simplified_mwall.rds"))
v <- c(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)
parcel_ids  <- 1:17
one_hot <- matrix(0, nrow = length(v), ncol = length(parcel_ids))
for (i in seq_along(parcel_ids)) {
  one_hot[, i] <- as.integer(v == parcel_ids[i])
}
mat <- dice_overlap(X=one_hot)
labs <- rownames(parcellation$meta$cifti$labels$parcels)[parcellation$meta$cifti$labels$parcels$Key > 0]
p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) +
  theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"), 
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))


# Make sure your tiles are square and labels readable
p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "Yeo17_parcellation_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)






# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "Yeo17_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "Yeo17_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6) 


# MSC
parcellation <- readRDS(file.path(dir_data, "outputs", "parcellations", "MSC_parcellation.rds"))
name = "MSC"
v <- c(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)
parcel_ids  <- 1:17
one_hot <- matrix(0, nrow = length(v), ncol = length(parcel_ids))

one_hot = one_hot[, order[[name]]$ic_order]

for (i in seq_along(parcel_ids)) {
  one_hot[, i] <- as.integer(v == parcel_ids[i])
}
mat <- dice_overlap(X=one_hot)
tab  <- parcellation$meta$cifti$labels$`Column number`
labs <- rownames(tab)[tab$Key > 0]
labs <- labs[order[[name]]$ic_order]

p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) +
  theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

 ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "MSC_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6) 

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "MSC_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none"
  )

fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),
  height = grid::unit(12, "cm")
)

gridExtra::grid.arrange(fixed_size_plot)

ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "MSC_parcellation_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)



# PROFUMO
parcellation <- readRDS(file.path(dir_data, "outputs", "parcellations", "PROFUMO_simplified_mwall.rds"))
name = "PROFUMO"
v <- rbind(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)           
sd <- apply(v, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(v), ncol(v), byrow = TRUE)
one_hot <- ifelse(abs(v) >= sd_mat, 1L, 0L)

one_hot = one_hot[, order[[name]]$ic_order]

mat <- dice_overlap(X=one_hot)
labs <- paste0("Network ", 1:12)
labs <- labs[order[[name]]$ic_order]
p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) + theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6)

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_parcellation_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)





# GICA15
parcellation <- read_cifti(file.path(dir_data, "inputs", "GICA15.dscalar.nii"))
name = "GICA15"
v <- rbind(
  parcellation$data$cortex_left,
  parcellation$data$cortex_right
)           
sd <- apply(v, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(v), ncol(v), byrow = TRUE)
one_hot <- ifelse(abs(v) >= sd_mat, 1L, 0L)
one_hot = one_hot[, order[[name]]$ic_order]
mat <- dice_overlap(X=one_hot)
labs <- paste0("IC", 1:15)
labs <- labs[order[[name]]$ic_order]
p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) + theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  )+scale_fill_gradientn(colours = viridisLite::mako(1024, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA15_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6) 

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA15_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )



p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "GICA15_parcellation_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)







# GICA25
# parcellation <- read_cifti(file.path(dir_data, "inputs", "GICA25.dscalar.nii"))
# v <- rbind(
#   parcellation$data$cortex_left,
#   parcellation$data$cortex_right
# )           
# sd <- apply(v, 2, sd, na.rm = TRUE)
# sd_mat  <- matrix(2 * sd, nrow(v), ncol(v), byrow = TRUE)
# one_hot <- ifelse(abs(v) >= sd_mat, 1L, 0L)
# mat <- dice_overlap(X=one_hot)
# labs <- paste0("IC", 1:25)
# p <- plot_FC_gg(
#   mat,
#   labs      = labs,
#   lim       = lim,
#   title="",
#   labs_margin_y = -15,
# ) + theme(
#     # legend.title = element_blank(),
#     # legend.text  = element_text(size = 14),
#     # legend.key.height = unit(2, "cm"),
#     # legend.key.width  = unit(0.6, "cm")
#     legend.position="none"
#   ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# # ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA25_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6)


# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA25_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


# # GICA50
# parcellation <- read_cifti(file.path(dir_data, "inputs", "GICA50.dscalar.nii"))
# v <- rbind(
#   parcellation$data$cortex_left,
#   parcellation$data$cortex_right
# )           
# sd <- apply(v, 2, sd, na.rm = TRUE)
# sd_mat  <- matrix(2 * sd, nrow(v), ncol(v), byrow = TRUE)
# one_hot <- ifelse(abs(v) >= sd_mat, 1L, 0L)
# mat <- dice_overlap(X=one_hot)
# labs <- paste0("IC", 1:50)
# p <- plot_FC_gg(
#   mat,
#   labs      = labs,
#   lim       = lim,
#   title="",
#   labs_margin_y = -8,
# ) + theme(
#     # legend.title = element_blank(),
#     # legend.text  = element_text(size = 14),
#     # legend.key.height = unit(2, "cm"),
#     # legend.key.width  = unit(0.6, "cm")
#     legend.position="none"
#   ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# # ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA50_parcellation_overlap.png"), plot = p, bg = "white", width=6, height=6)

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA50_parcellation_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )

########## PRIORS ##########


# Yeo17
prior <- readRDS(file.path(dir_project, "priors", "Yeo17", "prior_combined_Yeo17_noGSR.rds"))

# thresholded at 2SD
mean <- prior$prior$mean
sd  <- apply(mean, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
mask <- abs(mean) >= sd_mat
ones <- matrix(1L, nrow(mean), ncol(mean))
one_hot <- ones * mask

mat <- dice_overlap(X=one_hot)

labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]

p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) +
  theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "Yeo17_prior_overlap.png"), plot = p, bg = "white", width=6, height=6) 


# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "Yeo17_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "Yeo17_prior_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)




# MSC
prior <- readRDS(file.path(dir_project, "priors", "MSC", "prior_combined_MSC_noGSR.rds"))
name = "MSC"
mean <- prior$prior$mean[,2:18]
sd  <- apply(mean, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
mask <- abs(mean) >= sd_mat
ones <- matrix(1L, nrow(mean), ncol(mean))
one_hot <- ones * mask
one_hot = one_hot[, order[[name]]$ic_order]

mat <- dice_overlap(X=one_hot)
labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]
labs <- labs[order[[name]]$ic_order]

p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) +
  theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "MSC_prior_overlap.png"), plot = p, bg = "white", width=6, height=6) 

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "MSC_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "MSC_prior_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)





# PROFUMO
prior <- readRDS(file.path(dir_project, "priors", "PROFUMO", "prior_combined_PROFUMO_noGSR.rds"))
name = "PROFUMO"
mean <- prior$prior$mean
sd  <- apply(mean, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
mask <- abs(mean) >= sd_mat
ones <- matrix(1L, nrow(mean), ncol(mean))
one_hot <- ones * mask
one_hot = one_hot[, order[[name]]$ic_order]
mat <- dice_overlap(X=one_hot)
labs <- paste0("Network ", 1:ncol(prior$prior$mean))
labs <- labs[order[[name]]$ic_order]

p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) +
  theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_prior_overlap.png"), plot = p, bg = "white", width=6, height=6) 

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )


p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "PROFUMO_prior_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)


# GICA15 
prior <- readRDS(file.path(dir_project, "priors", "GICA15", "prior_combined_GICA15_noGSR.rds"))
name = "GICA15"
mean <- prior$prior$mean
sd  <- apply(mean, 2, sd, na.rm = TRUE)
sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
mask <- abs(mean) >= sd_mat
ones <- matrix(1L, nrow(mean), ncol(mean))
one_hot <- ones * mask
one_hot = one_hot[, order[[name]]$ic_order]
mat <- dice_overlap(X=one_hot)
labs <- paste0("IC", 1:15)
labs <- labs[order[[name]]$ic_order]
p <- plot_FC_gg(
  mat,
  labs      = labs,
  lim       = lim,
  title="",
  labs_margin_y = -14,
) + theme(
    # legend.title = element_blank(),
    # legend.text  = element_text(size = 14),
    # legend.key.height = unit(2, "cm"),
    # legend.key.width  = unit(0.6, "cm")
    legend.position="none"
  ) + scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA15_prior_overlap.png"), plot = p, bg = "white", width=6, height=6) 

# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))
# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA15_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )



p <- p + coord_fixed() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  )

# --- lock the actual matrix (panel) size ---
fixed_size_plot <- egg::set_panel_size(
  p = p,
  width  = grid::unit(12, "cm"),  # adjust to taste
  height = grid::unit(12, "cm")
)

# optional: preview
gridExtra::grid.arrange(fixed_size_plot)

# --- save ---
ggplot2::ggsave(
  filename = file.path(dir_data, "outputs", "dice_overlap", "GICA15_prior_overlap.png"),
  plot = fixed_size_plot,
  bg = "white",
  width = 8, height = 8, dpi = 300
)




# GICA25
# prior <- readRDS(file.path(dir_project, "priors", "GICA25", "prior_combined_GICA25_noGSR.rds"))
# mean <- prior$prior$mean
# sd  <- apply(mean, 2, sd, na.rm = TRUE)
# sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
# mask <- abs(mean) >= sd_mat
# ones <- matrix(1L, nrow(mean), ncol(mean))
# one_hot <- ones * mask
# mat <- dice_overlap(X=one_hot)
# labs <- paste0("IC", 1:25)
# p <- plot_FC_gg(
#   mat,
#   labs      = labs,
#   lim       = lim,
#   title="",
#   labs_margin_y = -15,
# ) + theme(
#     # legend.title = element_blank(),
#     # legend.text  = element_text(size = 14),
#     # legend.key.height = unit(2, "cm"),
#     # legend.key.width  = unit(0.6, "cm")
#     legend.position="none"
#   )+ scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# # ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA25_prior_overlap.png"), plot = p, bg = "white", width=6, height=6) 


# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))

# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA25_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )

# GICA50
# prior <- readRDS(file.path(dir_project, "priors", "GICA50", "prior_combined_GICA50_noGSR.rds"))
# mean <- prior$prior$mean
# sd  <- apply(mean, 2, sd, na.rm = TRUE)
# sd_mat  <- matrix(2 * sd, nrow(mean), ncol(mean), byrow = TRUE)
# mask <- abs(mean) >= sd_mat
# ones <- matrix(1L, nrow(mean), ncol(mean))
# one_hot <- ones * mask
# mat <- dice_overlap(X=one_hot)
# labs <- paste0("IC", 1:50)
# p <- plot_FC_gg(
#   mat,
#   labs      = labs,
#   lim       = lim,
#   title="",
#   labs_margin_y = -8,
# ) + theme(
#     # legend.title = element_blank(),
#     # legend.text  = element_text(size = 14),
#     # legend.key.height = unit(2, "cm"),
#     # legend.key.width  = unit(0.6, "cm")
#     legend.position="none"
#   )+ scale_fill_gradientn(colours = viridisLite::mako(256, direction = -1))

# # ggplot2::ggsave(file.path(dir_data, "outputs", "dice_overlap", "GICA50_prior_overlap.png"), plot = p, bg = "white", width=6, height=6)


# p_fixed <- set_panel_size(p, width = grid::unit(4, "in"), height = grid::unit(4, "in"))

# ggplot2::ggsave(
#   file.path(dir_data, "outputs", "dice_overlap", "GICA50_prior_overlap.png"),
#   plot = p_fixed,
#   bg = "white",
#   dpi = 300
# )

# #  + scale_fill_gradientn(colours = viridis::inferno(256, direction = -1))