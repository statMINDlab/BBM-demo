# Balance sex within age groups

# Load FD-filtered lists
valid_LR_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_LR_subjects_FD.csv"))$subject_id
valid_RL_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_RL_subjects_FD.csv"))$subject_id
valid_combined_subjects_FD <- read.csv(file.path(dir_data, "priors", "filtering", "valid_combined_subjects_FD.csv"))$subject_id

# Load unrelated-filtered lists if they exist
if (file.exists(file.path(dir_data, "priors", "filtering", "valid_LR_subjects_unrelated.csv"))) {
  valid_LR_subjects_unrelated <- read.csv(file.path(dir_data, "priors", "filtering", "valid_LR_subjects_unrelated.csv"))$subject_id
  valid_RL_subjects_unrelated <- read.csv(file.path(dir_data, "priors", "filtering", "valid_RL_subjects_unrelated.csv"))$subject_id
  valid_combined_subjects_unrelated <- read.csv(file.path(dir_data, "priors", "filtering", "valid_combined_subjects_unrelated.csv"))$subject_id
}

for (encoding in c("LR", "RL", "combined")) {
    
    if (file.exists(file.path(dir_data, "priors", "filtering", "valid_LR_subjects_unrelated.csv"))) {
        filtered_subjects <- HCP_unrestricted[HCP_unrestricted$Subject %in% get(sprintf("valid_%s_subjects_unrelated", encoding)), ]
    } else {
        # If skipping step 2 (no access to restricted data / not filtering by unrelated)
        filtered_subjects <- HCP_unrestricted[HCP_unrestricted$Subject %in% get(sprintf("valid_%s_subjects_FD", encoding)), ]
    }

    filtered_subjects$Age <- as.factor(filtered_subjects$Age)

    for (aa in seq_along(levels(filtered_subjects$Age))) {
        age_group <- levels(filtered_subjects$Age)[aa]
        subset_age <- subset(filtered_subjects, Age == age_group)
        gender_counts <- sort(table(subset_age$Gender))

        # Balance only if both genders exist
        if (length(gender_counts) == 2 && diff(gender_counts) != 0) {
            overrepresented_sex <- names(gender_counts)[2]
            to_sample_idx <- which(filtered_subjects$Age == age_group & filtered_subjects$Gender == overrepresented_sex)
            n_remove <- diff(gender_counts)
            subjects_to_remove_idx <- sample(to_sample_idx, n_remove)
            filtered_subjects <- filtered_subjects[-subjects_to_remove_idx, ]
        }
    }

    subject_ids_balanced <- filtered_subjects$Subject

    # Save new list
    if (encoding == "LR") {
        valid_LR_subjects_balanced <- subject_ids_balanced
    } else if (encoding == "RL") {
        valid_RL_subjects_balanced <- subject_ids_balanced
    } else {
        valid_combined_subjects_balanced <- subject_ids_balanced
    }  
}

write.csv(data.frame(subject_id=valid_LR_subjects_balanced), file = file.path(dir_data, "priors", "filtering", "valid_LR_subjects_balanced.csv"), row.names = FALSE)
write.csv(data.frame(subject_id=valid_RL_subjects_balanced), file = file.path(dir_data, "priors", "filtering", "valid_RL_subjects_balanced.csv"), row.names = FALSE)
write.csv(data.frame(subject_id=valid_combined_subjects_balanced), file = file.path(dir_data, "priors", "filtering", "valid_combined_subjects_balanced.csv"), row.names = FALSE)

saveRDS(valid_LR_subjects_balanced, file = file.path(dir_data, "priors", "filtering", "valid_LR_subjects_balanced.rds"))
saveRDS(valid_RL_subjects_balanced, file = file.path(dir_data, "priors", "filtering", "valid_RL_subjects_balanced.rds"))
saveRDS(valid_combined_subjects_balanced, file = file.path(dir_data, "priors", "filtering", "valid_combined_subjects_balanced.rds"))