
# # Plots Functional Connectivity (FC) priors for each prior using both the Cholesky and Inverse-Wishart parameterization
# #latest

# remove.packages("fMRItools")
# devtools::install_github("mandymejia/fMRItools", "7.0")
# library(fMRItools)

# remove.packages("BayesBrainMap")
# devtools::install_github("mandymejia/BayesBrainMap", "2.0")
# library(BayesBrainMap)

library(ggplot2)

# prior_files <- list.files(file.path(dir_project, "priors"), recursive = TRUE, full.names = TRUE)

# prior_files <- prior_files[grepl("combined", prior_files) & grepl("GSR", prior_files) & grepl("GICA15", prior_files)]

prior_files <- "~/Documents/GitHub/BBM-demo/data/priors/Yeo17_smooth-5/prior_combined_Yeo17_noGSR_local.rds"

get_prior_title <- function(base_name, encoding) {
  gsr <- if (grepl("noGSR", base_name)) "noGSR" else "GSR"

  if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
    return(paste0("Yeo17 ", gsr))
  } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
    return(paste0("MSC ", gsr))
  } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
    return(paste0("PROFUMO ", gsr))
  }

  ic_match <- regmatches(base_name, regexpr("GICA\\d+", base_name))
  nIC <- as.numeric(gsub("GICA", "", ic_match))

  paste0("GICA ", nIC, " ", gsr)
}

for (file in prior_files) {
    cat("Processing prior:", file, "\n")
    prior <- readRDS(file)
    
    base_name <- tools::file_path_sans_ext(basename(file))

    cat("Processing prior:", base_name, "\n")

    # LABELS
    if (grepl("Yeo17", base_name, ignore.case = TRUE)) {
      labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]
    } else if (grepl("MSC", base_name, ignore.case = TRUE)) {
      # change FC dim
      prior$prior$FC$Chol$mean <- prior$prior$FC$Chol$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$Chol$var  <- prior$prior$FC$Chol$var[2:18, 2:18, drop=FALSE]
      prior$prior$FC$IW$mean <- prior$prior$FC$IW$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$IW$var  <- prior$prior$FC$IW$var[2:18, 2:18, drop=FALSE]
      prior$prior$FC$empirical$mean <- prior$prior$FC$empirical$mean[2:18, 2:18, drop=FALSE]
      prior$prior$FC$empirical$var  <- prior$prior$FC$empirical$var[2:18, 2:18, drop=FALSE]
      labs <- rownames(prior$template_parc_table)[prior$template_parc_table$Key > 0]
    } else if (grepl("PROFUMO", base_name, ignore.case = TRUE)) {
      labs <- paste0("Network ", 1:12)
    } else {
      labs <- paste0("IC", 1:dim(prior$prior$mean)[2])
    }

    parts <- strsplit(base_name, "_")[[1]]
    encoding <- parts[2]      
    parcellation <- parts[3]   
    gsr_status <- parts[4]  

    out_dir <- file.path(dir_data, "priors", parcellation, "plots_FC")
    cat("out_dir:", out_dir, "\n")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Number of ICs
    Q <- dim(prior$prior$mean)[2]
    plot_title <- get_prior_title(base_name, encoding)
    p1 <- plot_FC_gg(
      prior$prior$FC$Chol$mean,
      labs      = labs,
      lim = c(-0.8, 0.8),
      title ="",
    #   labs_margin_y = -10
    ) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_Cholesky_mean.png")), plot = p1, bg = "white", width = 8, height = 8, dpi = 300)

    p2 <- plot_FC_gg(
      sqrt(prior$prior$FC$Chol$var),
      labs      = labs,
      lim = c(0, 0.4),
      title ="",
    #   labs_margin_y = -10
    )  +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_Cholesky_sd.png")), plot = p2, bg = "white", width = 8, height = 8, dpi = 300) 

    p3 <- plot_FC_gg(
      prior$prior$FC$IW$mean,
      labs      = labs,
      lim = c(-0.8, 0.8),
      title ="",
    #   labs_margin_y = -10
    )  +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_IW_mean.png")), plot = p3, bg = "white", width = 8, height = 8, dpi = 300)

    p4 <- plot_FC_gg(
      sqrt(prior$prior$FC$IW$var),
      labs      = labs,
      lim = c(0, 0.4),
      title ="",
    #   labs_margin_y = -10
    )  +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_IW_sd.png")), plot = p4, bg = "white", width = 8, height = 8, dpi = 300)

    p5 <- plot_FC_gg(
      prior$prior$FC$empirical$mean,
      labs      = labs,
      lim = c(-0.8, 0.8),
      title ="",
    #   labs_margin_y = -10
    )  +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_Empirical_mean.png")), plot = p5, bg = "white", width = 8, height = 8, dpi = 300)

    p6 <- plot_FC_gg(
      sqrt(prior$prior$FC$empirical$var),
      labs      = labs,
      lim = c(0, 0.4),
      title ="",
    #   labs_margin_y = -10
    )  +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 14),
    legend.key.height = unit(2, "cm"), 
    legend.key.width  = unit(0.6, "cm")
  )
    ggplot2::ggsave(file.path(out_dir, paste0(base_name, "_FC_Empirical_sd.png")), plot = p6, bg = "white", width = 8, height = 8, dpi = 300)
     

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

