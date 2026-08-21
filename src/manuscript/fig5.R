# Definitions
subid <- "100408"
sesid <- "REST2"
smoothing <- 5

# bMap
Q <- 17
# relevant component  
show_idx <- 14

# build path
dir_manuscript <- file.path(dir_project, "manuscript")

#bMap <- readRDS(file.path(dir_data, "outputs", "brain_map", "Yeo17", "four_sessions", "brainMap_4sessions.rds"))
bMap <- readRDS(file.path(dir_manuscript, "outputs", "brain_map", paste0("sub-", subid, "_ses-", sesid),
                          paste0("sub-", subid, "_ses-", sesid, "_brainmap_smoothed-", smoothing, "mm.rds")))

# Define prior
prior_fname <- if (!is.null(smoothing)) file.path(dir_data, paste0("priors/Yeo17_smooth-", smoothing, "/prior_combined_Yeo17_noGSR_local.rds")) else file.path(dir_data, paste0("priors/Yeo17/prior_combined_Yeo17_noGSR_local.rds"))
prior <- readRDS(prior_fname)
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)


label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == 14]
fname <- file.path(dir_manuscript, "Figure5", paste0("posterior_Yeo17_sub-", subid, "_ses-", sesid, "_", label_name))
dir.create(file.path(dir_manuscript, "Figure5"), recursive = TRUE)
# plot(bMap, idx = show_idx, stat = "mean", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, zlim=c(-round(3/Q, 3), round(3/Q, 3)))
# plot(bMap, idx = show_idx, stat = "se", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname,zlim=c(0, round(2/Q, 3)))


plot(bMap, idx = show_idx, stat = "mean", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, zlim=c(-0.4,0.4))
plot(bMap, idx = show_idx, stat = "se", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, zlim=c(0.021, 0.023))


# files_written <- plot(bMap, idx = show_idx, stat = "se", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, legend_fname = file.path("~/Desktop/Figure5/GICA_legend.png"), zlim=c(0, 0.040))
# 
# map_png    <- files_written[1]  
# legend_png <- files_written[2] 
# 
# ciftiTools::view_comp(
#     img    = map_png,
#     title=NULL,
#     legend = legend_png,
#     legend_height = 0.4,
#     fname  = file.path("~/Desktop", "Figure5", "Yeo17_bMap_Default_mean.png")
# )
# 
# fname <- file.path("~/Desktop/Figure5", paste0("posterior_Yeo17_", label_name))
# files_written <- plot(bMap, idx = show_idx, stat = "se", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, legend_fname = file.path("~/Desktop/Figure5/GICA_legend.png"),zlim=c(0, round(2/Q, 3)))
#  
# 
# map_png    <- files_written[1]  
# legend_png <- files_written[2] 
# 
# ciftiTools::view_comp(
#     img    = map_png,
#     title=NULL,
#     legend = legend_png,
#     legend_height = 0.4,
#     fname  = file.path("~/Desktop", "Figure5", "Yeo17_bMap_Default_se.png")
# )
# 
# bMap <- readRDS(file.path(dir_data, "outputs", "brain_map", "Yeo17", "four_sessions", "brainMap_4sessions.rds"))

eng <- engagements(bMap, z = c(0:3))
fname = file.path(dir_manuscript, "Figure5", paste0("posterior_Yeo17_sub-", subid, "_ses-", sesid, "_", label_name))
plot(eng, idx = show_idx, stat = "engaged", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname) 

# files_written <- plot(eng, idx = 14, stat = "engaged", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname, legend_fname = file.path("~/Desktop/Figure5/GICA_legend.png"), zlim=c(-0.2, 0.2)) 

# map_png    <- files_written[1]  
# legend_png <- files_written[2] 

# ciftiTools::view_comp(
#     img    = map_png,
#     title=NULL,
#     legend = legend_png,
#     legend_height = 0.5,
#     fname  = file.path("~/Desktop", "Figure5", "Yeo17_bMap_Default_eng.png")
# )







