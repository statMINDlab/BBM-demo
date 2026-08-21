# Yeo 17 -> idx 14
# GICA 15, 25 idx 2
# MSC -> 2 (Q=17)
# PROFUMO -> 8 (Q=12)

out_dir <- "~/Documents/GitHub/BBM-demo/manuscript"

# set var_method for exporting figures 
variance_method <- "non-negative"

# Select Parcellation
parcellation <- "Yeo17"
prior_file <- paste0("~/Documents/GitHub/BBM-demo/data/priors/", parcellation, "_smooth-5/prior_combined_", parcellation, "_noGSR.rds")
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 25
idx <- 2
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
    prior,
    what = "maps",
    stat  = "mean",
    idx   = idx,
    title = "",
    cex.title = 1e-6,
    legend_embed = FALSE,   
    zlim = c(-3/Q, 3/Q),                              
    fname = file.path(out_dir, "Figure2", "Yeo17_IC14.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "Yeo17_IC14.png"),
)


# GICA15
prior_file <- file.path(dir_data, "priors", sprintf("GICA%d_smooth-5/prior_combined_GICA%d_noGSR.rds", 15, 15))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 15
idx <- 2
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA15_IC2.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA15_IC2.png"),
)

# GICA25
prior_file <- file.path(dir_data, "priors", sprintf("GICA%d_smooth-5/prior_combined_GICA%d_noGSR.rds", 25, 25))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 25
idx <- 2
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA25_IC2.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA25_IC2.png"),
)

# MSC
prior_file <- file.path(dir_data, "priors", sprintf("MSC_smooth-5/prior_combined_MSC_noGSR_local.rds"))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 18
idx <- 2

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "MSC_IC2.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "MSC_IC2.png"),
)

# PROFUMO

prior_file <- file.path(dir_data, "priors", sprintf("PROFUMO_smooth-5/prior_combined_PROFUMO_noGSR_local.rds"))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 12
idx <- 8
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "PROFUMO_IC8.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "PROFUMO_IC8.png"),
)

profumo_mw <- readRDS(file.path(dir_data, "outputs", "parcellations", "PROFUMO_simplified_mwall.rds"))

plot(
        profumo_mw,
        idx   = 8,
        title = "",cex.title = 1e-6,
        legend_embed = FALSE,                                    
        fname = file.path(out_dir, "Figure2", "PROFUMO_8.png")
)





files_written <- plot(
        profumo_mw,
        idx   = 8,
        title = "",cex.title = 1e-6,
        legend_embed = FALSE,                                 
        fname = file.path(out_dir, "PROFUMO"),      
        legend_fname = file.path(out_dir, "PROFUMO_legend.png")
)

    # map_png    <- files_written$surface[1]  
    # legend_png <- files_written$surface[2] 


map_png    <- files_written[1]  
legend_png <- files_written[2] 

ciftiTools::view_comp(
    img    = map_png,title=NULL,
    legend = legend_png,
    legend_height = 0.5,
    fname  = file.path(out_dir, "Figure2", "PROFUMO_8.png")
) 












# Set the prior parameters
Q <- 25
idx <- 2
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA25_IC2.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "GICA25_IC2.png"),
)

# MSC
prior_file <- file.path(dir_data, "priors", sprintf("MSC_smooth-5/prior_combined_MSC_noGSR_local.rds"))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 18
idx <- 14

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "MSC_IC2.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "MSC_IC2.png"),
)

# PROFUMO

prior_file <- file.path(dir_data, "priors", sprintf("PROFUMO/prior_combined_PROFUMO_noGSR_local.rds"))
prior <- readRDS(prior_file)

# Set the prior parameters
Q <- 12
idx <- 8
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)

plot(
  prior,
  what = "maps",
  stat  = "mean",
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(-3/Q, 3/Q),                              
  fname = file.path(out_dir, "Figure2", "PROFUMO_IC8.png"),
)

plot(
  prior,
  what = "maps",
  stat  = "sd", var_method = variance_method,
  idx   = idx,
  title = "",
  cex.title = 1e-6,
  legend_embed = FALSE,   
  zlim = c(0, 2/Q),                              
  fname = file.path(out_dir, "Figure2", "PROFUMO_IC8.png"),
)

profumo_mw <- readRDS(file.path(dir_data, "outputs", "parcellations", "PROFUMO_simplified_mwall.rds"))

plot(
        profumo_mw,
        idx   = 8,
        title = "",cex.title = 1e-6,
        legend_embed = FALSE,                                    
        fname = file.path(out_dir, "Figure2", "PROFUMO_8.png")
)





files_written <- plot(
        profumo_mw,
        idx   = 8,
        title = "",cex.title = 1e-6,
        legend_embed = FALSE,                                 
        fname = file.path(out_dir, "PROFUMO"),      
        legend_fname = file.path(out_dir, "PROFUMO_legend.png")
)

    # map_png    <- files_written$surface[1]  
    # legend_png <- files_written$surface[2] 


map_png    <- files_written[1]  
legend_png <- files_written[2] 

ciftiTools::view_comp(
    img    = map_png,title=NULL,
    legend = legend_png,
    legend_height = 0.5,
    fname  = file.path(out_dir, "Figure2", "PROFUMO_8.png")
) 













