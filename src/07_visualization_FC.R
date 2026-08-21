
# # Plots Functional Connectivity (FC) priors for each prior using both the Cholesky and Inverse-Wishart parameterization
# #latest

# remove.packages("fMRItools")
# devtools::install_github("mandymejia/fMRItools", "8.0") # Version required for match_nets()
library(fMRItools)

# remove.packages("BayesBrainMap")
# devtools::install_github("mandymejia/BayesBrainMap", "2.0")
# library(BayesBrainMap)

library(ggplot2)

prior_files <- list.files(file.path(dir_data, "priors"), recursive = TRUE, full.names = TRUE, pattern = "*GSR_local.rds")

source(file.path(dir_project, "src", "06_best_match_IC.R"))


get_prior_title <- function(base_name, i, prior, encoding, gsr_status) {
  
  gsr <- if (grepl("noGSR", base_name)) "noGSR" else "GSR"
  
  if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
    return(paste0("Yeo17 ", gsr))
  } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
    return(paste0("MSC ", gsr))
  } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
    return(paste0("PROFUMO ", gsr))
  } else if (grepl("NMF", base_name, ignore.case = TRUE)) {
    return(paste0("NMF ", gsr))
  } 
  ic_match <- regmatches(base_name, regexpr("GICA\\d+", base_name))
  
  nIC <- as.numeric(gsub("GICA", "", ic_match))
  title_str <- paste0("GICA ", nIC, " - Component ", i)
  
  return(title_str)
}

plot_fc_all <- function(template, encoding, labs, out_dir, plot_title, base_name) {
  FC <- template$prior$FC

  legend_theme <- theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )

  plots <- list(
    list(mat = FC$Chol$mean,            lim = c(-0.8, 0.8), suffix = "Cholesky_mean",  title = paste0(plot_title, " Cholesky FC Prior Mean")),
    list(mat = sqrt(FC$Chol$var),       lim = c(0,    0.4), suffix = "Cholesky_sd",    title = paste0(plot_title, " Cholesky FC Prior SD")),
    list(mat = FC$IW$mean,              lim = c(-0.8, 0.8), suffix = "IW_mean",        title = paste0(plot_title, " Inverse-Wishart FC Prior Mean")),
    list(mat = sqrt(FC$IW$var),         lim = c(0,    0.4), suffix = "IW_sd",          title = paste0(plot_title, " Inverse-Wishart FC Prior SD")),
    list(mat = FC$empirical$mean,       lim = c(-0.8, 0.8), suffix = "Empirical_mean", title = paste0(plot_title, " Empirical FC Prior Mean")),
    list(mat = sqrt(FC$empirical$var),  lim = c(0,    0.4), suffix = "Empirical_sd",   title = paste0(plot_title, " Empirical FC Prior SD"))
  )

  for (p in plots) {
    plt <- plot_FC_gg(
      p$mat,
      labs          = labs,
      lim           = p$lim,
      labs_margin_y = -10,
      title         = p$title
    ) + legend_theme
    ggplot2::ggsave(
      file.path(out_dir, paste0(base_name, "_FC_", p$suffix, ".png")),
      plot = plt, bg = "white"
    )
  }
}

for (file in prior_files) {
    cat("Processing prior:", file, "\n")
    prior <- readRDS(file)

    base_name <- tools::file_path_sans_ext(basename(file))

    cat("Processing prior:", base_name, "\n")

    # LABELS
    if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
      labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]
      name = "Yeo17"
    } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
      # change FC dim
      prior$prior$FC$Chol$mean       <- prior$prior$FC$Chol$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$Chol$var        <- prior$prior$FC$Chol$var[2:18, 2:18, drop=FALSE]
      prior$prior$FC$IW$mean         <- prior$prior$FC$IW$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$IW$var          <- prior$prior$FC$IW$var[2:18, 2:18, drop=FALSE]
      prior$prior$FC$empirical$mean  <- prior$prior$FC$empirical$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$empirical$var   <- prior$prior$FC$empirical$var[2:18, 2:18, drop=FALSE]
      labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]
      name = "MSC"
    } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
      labs <- paste0("Network ", 1:12)
      name = "PROFUMO"
    } else if (grepl("NMF", base_name, ignore.case = TRUE)) {
      labs <- paste0("NMF ", pnc_labels)
      name = "NMF"
    } else {
      labs <- paste0("IC", 1:dim(prior$prior$mean)[2])
    }

    if (grepl("GICA15", base_name, ignore.case = TRUE)) {
      name = "GICA15"
    } else if (grepl("GICA25", base_name, ignore.case = TRUE)) {
      name = "GICA25"
    }

    parts <- strsplit(base_name, "_")[[1]]
    encoding    <- parts[2]
    parcellation <- parts[3]
    gsr_status  <- parts[4]

    
    out_dir <- file.path(dir_data, "priors", parcellation, "plots_FC")
    cat("out_dir:", out_dir, "\n")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    plot_title <- get_prior_title(base_name, encoding)

    ###############################################################
    # Reorder ICs to match Yeo's canonical order (from 06_best_match_IC.R)
    ###############################################################
    labs <- labs[order[[name]]$ic_order]
    prior$prior$FC$Chol$mean      <- prior$prior$FC$Chol$mean[order[[name]]$ic_order, order[[name]]$ic_order]
    prior$prior$FC$IW$mean        <- prior$prior$FC$IW$mean[order[[name]]$ic_order, order[[name]]$ic_order]
    prior$prior$FC$empirical$mean <- prior$prior$FC$empirical$mean[order[[name]]$ic_order, order[[name]]$ic_order]
    prior$prior$FC$Chol$var       <- prior$prior$FC$Chol$var[order[[name]]$ic_order, order[[name]]$ic_order]
    prior$prior$FC$IW$var         <- prior$prior$FC$IW$var[order[[name]]$ic_order, order[[name]]$ic_order]
    prior$prior$FC$empirical$var  <- prior$prior$FC$empirical$var[order[[name]]$ic_order, order[[name]]$ic_order]

    plot_fc_all(prior, encoding, labs, out_dir, plot_title, base_name)
}




# # ALL TOGETHER
# library(gridExtra)
# library(grid)
# plots <- list(
#   Chol_Mean = p1,
#   Chol_SD   = p2,
#   IW_Mean   = p3,
#   IW_SD     = p4,
#   Emp_Mean  = p5,
#   Emp_SD    = p6
# )
# combined <- grid.arrange(
#   plots$Chol_Mean, plots$Chol_SD, plots$IW_Mean,
#   plots$IW_SD,    plots$Emp_Mean, plots$Emp_SD,
#   nrow = 2, ncol = 3
# )

# ggsave(
#   file.path("~/Desktop", paste0(base_name, "_FC_ALL.png")),
#   combined,
#   width = 16, height = 9, dpi = 300, bg = "white"
# )
