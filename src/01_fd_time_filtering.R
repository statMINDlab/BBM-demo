# Filter Subjects by Sufficient fMRI Scan Duration

# Declare empty lists 
valid_LR_subjects_FD <- c()
valid_RL_subjects_FD <- c()
fd_flags <- list()

# Initialize table
fd_summary <- data.frame(
  subject = character(),
  session = character(),
  encoding = character(),
  mean_fd = numeric(),
  valid_time_sec = numeric()
)

# subject ids will be 

for (subject in subject_ids) {
    for (encoding in c("LR", "RL")) {

        session_pass <- c()

        for (session in c("REST1", "REST2")) {
            
            path <- sprintf("%s/%s/MNINonLinear/Results/rfMRI_%s_%s/Movement_Regressors.txt", dir_HCP, subject, session, encoding)

            # If file does not exist automatically false
            if (!file.exists(path)) {
                session_pass <- c(session_pass, FALSE)
                next
            }

            X <- as.matrix(read.table(path))

            fd <- FD(X=X[,1:6], lag=fd_lag_HCP, cutoff=fd_cutoff, TR_for_resp_filt=TR_HCP)

            mean_fd <- mean(fd$measure, na.rm = TRUE)
            
            # Use logical array to determine the valid volumes (below the cutoff, 0's in array)
            valid_volumes <- sum(!fd$outlier_flag)
            total_time_sec <- TR_HCP * valid_volumes

            # Check condition for filtering
            # No mean fd condition, only total time
            passed <- total_time_sec >= min_total_sec
            session_pass <- c(session_pass, passed)
            
            fd_summary <- rbind(fd_summary, data.frame(
                subject = subject,
                session = session,
                encoding = encoding,
                mean_fd = mean_fd,
                valid_time_sec = total_time_sec
            ))
            
            # for the list, subject id needs to be string
            subject_str <- as.character(subject)
            
            # Save fd object for future reference
            fd_flags[[encoding]][[subject_str]][[session]] <- fd$outlier_flag
        }
        # If both rest 1 and rest 2 meet both requirements, add to the list of valid for LR or RL
        if (all(session_pass)) {
            if (encoding == "LR") {
                valid_LR_subjects_FD <- c(valid_LR_subjects_FD, subject)
            } else if (encoding == "RL") {
                valid_RL_subjects_FD <- c(valid_RL_subjects_FD, subject)
            }
        }
    }
}

# Which subjects have valid LR and RL data?
valid_combined_subjects_FD <- intersect(valid_LR_subjects_FD, valid_RL_subjects_FD)

dir.create(file.path(dir_data, "outputs", "filtering"), recursive = TRUE, showWarnings = FALSE)

write.csv(data.frame(subject_id=valid_LR_subjects_FD), file = file.path(dir_data, "priors", "filtering", "valid_LR_subjects_FD.csv"), row.names = FALSE)
write.csv(data.frame(subject_id=valid_RL_subjects_FD), file = file.path(dir_data, "priors", "filtering", "valid_RL_subjects_FD.csv"), row.names = FALSE)
write.csv(data.frame(subject_id=valid_combined_subjects_FD), file = file.path(dir_data, "priors", "filtering", "valid_combined_subjects_FD.csv"), row.names = FALSE)
write.csv(fd_summary, file = file.path(dir_data, "priors", "filtering", "fd_summary.csv"), row.names = TRUE)

saveRDS(fd_flags, file.path(dir_data, "priors", "filtering", "fd_flags.rds"))


