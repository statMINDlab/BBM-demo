# Filter unrelated subjects

if (!file.exists(HCP_restricted_fname)) {
  warning("Skipping unrelated subject filtering: 'restricted_HCP.csv' not found in data_OSF/inputs")
} else {
  HCP_restricted <- read.csv(HCP_restricted_fname)

  valid_LR_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_LR_subjects_FD.csv"))$subject_id
  valid_RL_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_RL_subjects_FD.csv"))$subject_id
  valid_combined_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_combined_subjects_FD.csv"))$subject_id

  for (encoding in c("LR", "RL", "combined")) {
    # From restricted data get only the rows that correspond to the valid subject IDs after FD correction
    filtered_restricted <- HCP_restricted[HCP_restricted$Subject %in% get(sprintf("valid_%s_subjects_FD", encoding)), ]
    
    # Keep only one subject per Family ID
    filtered_unrelated <- filtered_restricted[!duplicated(filtered_restricted$Family_ID), ]

    # Get only subjects id from filtered_unrelated
    subject_ids_unrelated <- filtered_unrelated$Subject

    # Save new list
    if (encoding == "LR") {
      valid_LR_subjects_unrelated <- subject_ids_unrelated
    } else if (encoding == "RL") {
      valid_RL_subjects_unrelated <- subject_ids_unrelated
    } else {
      valid_combined_subjects_unrelated <- subject_ids_unrelated
    }
  }

  # Save intermediate data in personal directory (due to restricted data)
  write.csv(data.frame(subject_id=valid_LR_subjects_unrelated), file = file.path(dir_data, "priors", "filtering", "valid_LR_subjects_unrelated.csv"), row.names = FALSE)
  write.csv(data.frame(subject_id=valid_RL_subjects_unrelated), file = file.path(dir_data, "priors", "filtering", "valid_RL_subjects_unrelated.csv"), row.names = FALSE)
  write.csv(data.frame(subject_id=valid_combined_subjects_unrelated), file = file.path(dir_data, "priors", "filtering", "valid_combined_subjects_unrelated.csv"), row.names = FALSE)
}