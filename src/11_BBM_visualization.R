#manuscript_brainmap_visualization.R


# Initialize dependencies 
sourcedir = "~/Documents/GitHub/BayesianBrainMapping-priors/src"

# Setup up dependencies and parameters
source(file.path(sourcedir, "0_setup.R"))

get_prior_title <- function(base_name, i, prior, encoding, gsr_status) {
  
  if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
    label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]
    return(paste0("Yeo 17 Network ", label_name, " (#", i, ")"))
  } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
    label_name <- rownames(prior$template_parc_table)[i]
    return(paste0("MSC Network ", label_name, " (#", i-1, ")"))
  } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
    return(paste0("PROFUMO Network # ", i))
  }
  ic_match <- regmatches(base_name, regexpr("GICA\\d+", base_name))
  
  nIC <- as.numeric(gsub("GICA", "", ic_match))
  title_str <- paste0("GICA ", nIC, " - Component ", i)
  
  return(title_str)
}

################################### Set parameters to look-up RDS names. ###################################################################

manuscript_output_dir <- "~/Documents/GitHub/BayesianBrainMapping-priors/manuscript"
output_dir <- file.path(manuscript_output_dir, "outputs", "brain_map")
# create output directory if it does not exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set "fake" transparent colours for visualization.
transparent_colors <- c("#ff8080", "#8080ff", "#CC80B3")

# set subject and session
subject_ids <- c("100307", "100206", "100408", "100610", "101107", "103111") # example subject 
session_id <- c("REST1", "REST2")

# set parameters
encoding <- c("LR", "RL") 
smoothing <- 4 # in mm FWHM
scrubbing <- TRUE
GSR = FALSE
# Define prior path based on selected nIC
nIC <- brainMap_prior
prior_path <- if (nIC == 0) {
  file.path(dir_project, "priors", "Yeo17", "prior_combined_Yeo17_noGSR.rds")
} else if (nIC == 1) {
  file.path(dir_project, "priors", "MSC", "prior_combined_MSC_noGSR.rds")
} else if (nIC == 2) {
  file.path(dir_project, "priors", "PROFUMO", "prior_combined_PROFUMO_noGSR.rds")
} else {
  file.path(dir_project, "priors", sprintf("GICA%d", nIC), paste0("prior_combined_", sprintf("GICA%d", nIC), "_noGSR.rds"))
}

# get label names
prior <- readRDS(prior_path)
prior_name = "Yeo17"
prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)
Q <- dim(prior$prior$mean)[2]

############################## Start plotting #########################################################################################

for (subid in subject_ids){

  subject_engagements <- list()

  for (sesid in session_id) {

    # Get name

    cat("Subject: ", subid, " Session:", sesid, "\n")
    
    # Define base name for outputs
    base_name <- paste0("sub-", subid, "ses-", sesid, "_brainmap")
    
    # Define brain map output directory for the subject
    bm_dir <- file.path(output_dir, paste0("sub-", subid, "_ses-", sesid))
    
    if (scrubbing) {
      
      # check if scrubbing results already exist
      scrubbing_file <- file.path(bm_dir, paste0(base_name, "_scrubbing_results.rds"))
      }
    
    if (smoothing) {
      # add _smoothed-XXmm to the file names
      base_name <- paste0(base_name, "_smoothed-", smoothing, "mm")
    }

    # Read in RDS file  
    bMap = readRDS(file.path(bm_dir, paste0(base_name, ".rds")))

    # Calculate engagements for all networks
    z = c(1, 2, 3)
    eng <- id_engagements(
      bMap,
      z = z,
      method_p = "bonferroni"
    )
    
    # Plot brainmap and engagement for each network separately
    for (i in 1:Q){
      
        # Get label for component
        label_name = rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]          
        
        # plot brainMap scalar maps
        fname = file.path(dir_data, "../manuscript/plots", paste0("sub-", subid, "_ses-", sesid, "_", prior_name, "-", label_name, ".png"))
        plot(bMap, idx = i, stat = "mean", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname) 
    
        # Generate engagement map ############## FIGURE FOCAL ENGAGEMENT MAP ##############
        plot(eng, idx = i, stat = "engaged", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname) 


    }
    
    # Save engagement map at Z>1 for comparison between sessions

  subject_engagements[[sesid]] <- engagements(bMap, z = 1, method_p = "bonferroni")
    
  }

  # Make comparison between sessions ############## FIGURE TEST-RESTEST RELIABILITY ######

  # open both sessions of the subject
  # make cifti file with both sessions 
  comparison_cifti = (subject_engagements[[session_id[2]]]$engaged * 2 + subject_engagements[[session_id[1]]]$engaged)
  # Medialwall is expected as -1
  comparison_cifti$data$cortex_left[comparison_cifti$data$cortex_left <= 0] = NA
  comparison_cifti$data$cortex_right[comparison_cifti$data$cortex_right <= 0] = NA
  
  for (i in 1:Q){
    
    # Get label for component
    label_name = rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]          
    
    # plot brainMap scalar maps
    fname = file.path(dir_data, "../manuscript/plots", paste0("sub-", subid, "_ses-", sesid, "_", prior_name, "-", label_name))
    plot(comparison_cifti, idx = i, alpha = 1, colors = transparent_colors, color_mode = "qualitative", bg = "white", NA_color="white", fname=fname)

    # plot comparison
    #plot(comparison_cifti, idx = 14, stat = "engaged", title = "", cex.title = 1e-6, legend_embed = FALSE, fname=fname)

  }


}
