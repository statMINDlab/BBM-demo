# Plots both the mean and standard deviation components for all priors

prior_files <- list.files(file.path(dir_data, "priors"), recursive = TRUE, full.names = TRUE)

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

for (file in prior_files) {
  prior <- readRDS(file)

  base_name <- tools::file_path_sans_ext(basename(file))

  parts <- strsplit(base_name, "_")[[1]]
  encoding <- parts[2]      
  parcellation <- parts[3]   
  gsr_status <- parts[4]  

  # If Yeo17 template, template_parc_table needs to be updated to only reflect the correct number of labels (17)
  if (grepl("Yeo17", base_name)) {
    prior$template_parc_table <- subset(prior$template_parc_table, prior$template_parc_table$Key > 0)
  }

  Q <- dim(prior$prior$mean)[2]
  # Save 4 images for each IC (cortical sd and mean, and subcortical sd and mean)
  for (i in 1:Q) {
    if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
      label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == i]  
      fname <- file.path(dir_data, "outputs", "priors_plots", parcellation, encoding, gsr_status, paste0(base_name, "_", label_name))
    } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
      label_name <- rownames(prior$template_parc_table)[i]
      fname <- file.path(dir_data, "outputs", "priors_plots", parcellation, encoding, gsr_status,
                         paste0(base_name, "_", label_name))
    } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
      # label_name <- rownames(prior$template_parc_table)[prior$template_parc_table$Key == i-1]  
      fname <- file.path(dir_data, "outputs", "priors_plots", parcellation, encoding, gsr_status, paste0(base_name, "_", i))
      print(fname)
    } else {
      fname <- file.path(dir_data, "outputs", "priors_plots",  parcellation, encoding, gsr_status, paste0(base_name, "_IC", i))
    }

    outdir <- dirname(fname)
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    title <- get_prior_title(base_name, i, prior, encoding, gsr_status)
    cat("Plotting prior ", i, "\n")

    # Plot mean
    plot(
      prior,
      stat = "mean",
      fname = fname,
      idx = i,
      title = title
    )

    # Plot standard deviation
    plot(
      prior,
      stat = "sd",
      fname = fname,
      idx = i,
      title = title
    )
  }
}
