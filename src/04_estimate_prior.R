# Estimate Priors using `estimate_prior()`

# Example func call: estimate_and_export_prior("LR", 15, FALSE, dir_data, TR_HCP)
# encoding is "LR" / "RL" / "combined"
# nIC is 15 / 25 / 50, 0 meaning it is going to use the Yeo17 parcellation, 1 meaning MSC parcellation, or 2 meaning PROFUMO parcellation
# GSR is TRUE / FALSE
estimate_and_export_prior <- function(
  encoding,
  nIC,
  GSR,
  dir_data,
  TR_HCP,
  usePar
) {
    # Get final list of subjects 
    final_subject_ids <- readRDS(file.path(dir_data, "priors", "filtering", sprintf("valid_%s_subjects_balanced.rds", encoding)))

    # Construct file paths
    if (encoding == "LR" | encoding == "RL") {
        BOLD_paths1 <- file.path(dir_HCP, 
                                final_subject_ids, 
                                sprintf("MNINonLinear/Results/rfMRI_REST1_%s/rfMRI_REST1_%s_Atlas_MSMAll_hp2000_clean.dtseries.nii", encoding, encoding))
        encoding1 = encoding
        session1 = "REST1"
    
        BOLD_paths2 <- file.path(dir_HCP, 
                                final_subject_ids, 
                                sprintf("MNINonLinear/Results/rfMRI_REST2_%s/rfMRI_REST2_%s_Atlas_MSMAll_hp2000_clean.dtseries.nii", encoding, encoding))
        encoding2 = encoding
        session2 = "REST2"
    } else {
        BOLD_paths1 <- file.path(dir_HCP, 
                                final_subject_ids, 
                                sprintf("MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii"))
        encoding1 = "LR"
        session1 = "REST1"
    
        BOLD_paths2 <- file.path(dir_HCP, 
                                final_subject_ids, 
                                sprintf("MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii"))
        encoding2 = "RL"
        session2 = "REST1"
    }

    parcellation <- if (nIC == 0) {
        "Yeo17"
    } else if (nIC == 1) {
        "MSC"
    } else if (nIC == 2) {
        "PROFUMO"
    } else {
        sprintf("GICA%d", nIC)
    }
    
    gsr_label <- ifelse(GSR, "GSR", "noGSR")
    save_dir <- file.path(dir_data, "priors", parcellation)
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

    cat(sprintf("Estimating prior for encoding: %s , parcellation: %s , %s, Using %s threads\n",encoding, parcellation, gsr_label, as.character(usePar)))
    
    # Start scrubbing procedure, keeping ten minutes below FD threshold.
    
    # get fd flags from inputs
    fd_flags = readRDS(file.path(dir_data, "priors", "filtering", "fd_flags.rds"))
    
    # Define number of volumes to keep
    keep_volumes <-floor(min_total_sec / TR_HCP)
    
    # make nested list fd_flags into tibble, for easy vectorization
    fd_tbl <- fd_flags %>%
      imap_dfr(function(enc_list, encoding) {
        enc_list %>%
          imap_dfr(function(subj_list, subject_id) {
            subj_list %>%
              imap_dfr(function(fd_vec, session) {
                tibble(
                  encoding = encoding,
                  subject  = subject_id,
                  session  = session,
                  fd       = list(fd_vec)
                )
              })
          })
      })
    
    # filtered tibble for final subject list 
    fd_tbl <- fd_tbl %>%
      filter(subject %in% final_subject_ids)
    
    # Obtain the indices to keep, all subjects represented by 10 minutes
    fd_tbl <- fd_tbl %>%
      mutate(
        fd_scrubbed = map(
          fd,
          ~ {
            over_threshold <- !.x
            (cumsum(over_threshold) <= keep_volumes) & over_threshold
          }
        )
      )
    
    # make sure that I am keeping the required amount of volumes for each subject
    stopifnot(all(unlist(lapply(fd_tbl$fd_scrubbed, sum)) == keep_volumes))
    
    # format scrub indices in BBM-friendly way
    scrub_BOLD1 <- fd_tbl %>%
      filter(encoding == encoding1 & session == session1) %>%
      select(fd_scrubbed)
    scrub_BOLD2 <- fd_tbl %>%
      filter(encoding == encoding2 & session == session2) %>%
      select(fd_scrubbed)
    scrub <- list(unlist(scrub_BOLD1, recursive = FALSE), 
                  unlist(scrub_BOLD2, recursive = FALSE))


    # FD scrubbing
    

    
    # Yeo17 parcellation
    if (nIC == 0) {

        GICA <- readRDS(file.path(dir_data, "templates", "Yeo17_simplified_mwall.rds"))

        # Include certain ICs (1:17 not 0 or -1 -> medial wall)
        valid_keys <- GICA$meta$cifti$labels[[1]]$Key
        inds <- valid_keys[valid_keys > 0]

        # Diego Debug print libpaths

        prior <- estimate_prior(
                BOLD = BOLD_paths1,
                BOLD2 = BOLD_paths2,
                template = GICA,
                GSR = GSR,
                TR = TR_HCP,
                hpf = 0.01,
                Q2 = 0,
                Q2_max = NULL,
                verbose = TRUE,
                inds = inds,
                brainstructures = c("left", "right"),
                drop_first = 15,
                scrub = scrub,
                usePar=usePar,
                wb_path = wb_path
            )
        
        # Save file
        saveRDS(prior, file.path(save_dir, sprintf("prior_%s_%s_%s.rds", encoding, parcellation, gsr_label)))

    # MSC
    } else if (nIC == 1) {

        GICA <- readRDS(file.path(dir_data, "templates", "MSC_parcellation.rds"))

        prior <- estimate_prior(
                BOLD = BOLD_paths1,
                BOLD2 = BOLD_paths2,
                template = GICA,
                GSR = GSR,
                TR = TR_HCP,
                hpf = 0.01,
                Q2 = 0,
                Q2_max = NULL,
                verbose = TRUE,
                brainstructures = c("left", "right"),
                drop_first = 15,
                scrub = scrub,
                usePar=usePar,
                wb_path = wb_path
            )

        # Save file
        saveRDS(prior, file.path(save_dir, sprintf("prior_%s_%s_%s.rds", encoding, parcellation, gsr_label)))
    
    # PROFUMO
    } else if (nIC == 2) {

        PROFUMO <- readRDS(file.path(dir_data, "templates", "PROFUMO_simplified_mwall.rds"))

        prior <- estimate_prior(
                BOLD = BOLD_paths1,
                BOLD2 = BOLD_paths2,
                template = PROFUMO,
                GSR = GSR,
                TR = TR_HCP,
                hpf = 0.01,
                Q2 = 0,
                Q2_max = NULL,
                verbose = TRUE,
                brainstructures = c("left", "right"),
                drop_first = 15,
                scrub = scrub,
                usePar=usePar,
                wb_path = wb_path
            )

        # Save file
        save_dir <- file.path(dir_data, "priors", "PROFUMO")
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
        saveRDS(prior, file.path(save_dir, sprintf("prior_%s_PROFUMO_%s.rds", encoding, gsr_label)))

    # GICA
    } else {

        GICA <- file.path(dir_data, "templates", sprintf("GICA%d.dscalar.nii", nIC))

        prior <- estimate_prior(
                BOLD = BOLD_paths1,
                BOLD2 = BOLD_paths2,
                template = GICA,
                GSR = GSR,
                TR = TR_HCP,
                hpf = 0.01,
                Q2 = 0,
                Q2_max = NULL,
                verbose = TRUE,
                drop_first = 15,
                scrub = scrub,
                usePar=usePar,
                wb_path = wb_path
                )

        # Save file
        saveRDS(prior, file.path(save_dir, sprintf("prior_%s_%s_%s.rds", encoding, parcellation, gsr_label)))
    }

    cat(sprintf("Saved prior for encoding: %s , parcellation: %s , %s\n",encoding, parcellation, gsr_label))
}
