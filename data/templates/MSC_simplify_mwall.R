out_dir <- file.path(dir_data, "outputs", "parcellations_plots", "MSC")
nNet <- 17
network_df <- data.frame(
  name = c(
    "Medial Wall",
    "Default", "Lat Vis", "Fronto-Par",
    "Med Vis", "Dors Attn", "Premotor",
    "Vent Attn", "Salience", "Cing-Operc",
    "Hand SM", "Face SM", "Auditory",
    "Ant MTL", "Post MTL", "Par Memory",
    "Context", "Foot SM"
  ),
  name_merged = c(
    "NA",
    "Default", "Visual", "Fronto-Par",
    "Visual", "Dors Attn", "Premotor",
    "Vent Attn", "Salience", "Cing-Operc",
    "Motor", "Motor", "Auditory",
    "Ant MTL", "Post MTL", "Par Memory",
    "Context", "Motor"
  ),
  name_merged_unique = c(
    "NA",
    "Default", "Vis", "Fronto-Par",
    "Vis2", "Dors Attn", "Premotor",
    "Vent Attn", "Salience", "Cing-Operc",
    "Motor", "Motor2", "Auditory",
    "Ant MTL", "Post MTL", "Par Memory",
    "Context", "Motor3"
  ),
  color = c(
    "#cccccc",
    "#ff0002", "#002bd0", "#fffb00",
    "#ffcf9a", "#1efc0a", "#feb1f5",
    "#48b9d8", "#000000", "#5b1689",
    "#0fffec", "#ffbb00", "#c17cc2",
    "#577fba", "#9dfc9a", "#0359ff",
    "#7f7f7f", "#548236"
  ),
  value = seq(0, nNet)
)

mparc <- read_cifti(file.path(dir_data, "inputs", "Networks_template.dscalar.nii"))
# mparc <- move_from_mwall(mparc, 0)

invalid_vals <- setdiff(unique(mparc$data$subcort), network_df$value)
mparc$data$subcort[mparc$data$subcort %in% invalid_vals] <- NA

levels_old <- network_df$value
labels <- network_df$name
colors <- network_df$color

mparc <- convert_xifti(
  mparc, to = "dlabel",
  levels_old = levels_old,
  levels = levels_old,
  labels = labels,
  colors = colors,
  add_white = FALSE
)

mparc <- remove_xifti(mparc, remove = "subcort")

plot(mparc,
    fname = file.path(out_dir, "MSC.png"),
    title = "MSC Network Parcellation"
)

# Save parcellation
saveRDS(
  mparc,
  file.path(dir_data, "outputs", "MSC_parcellation.rds")
)

# Plot every network separately

for (parc in 1:17) {
    
  parc_copy <- mparc
  curr_color <- network_df[network_df$value == parc, "color"]
  rgb_vals <- col2rgb(curr_color)
  labels_table <- parc_copy$meta$cifti$labels[[1]]
  labels_table[, 2:4] <- matrix(1, nrow = nrow(labels_table), ncol = 3)
  labels_table[labels_table$Key == parc, 2:4] <- t(rgb_vals / 255)
  parc_copy$meta$cifti$labels[[1]] <- labels_table
  label_name <- network_df[network_df$value == parc, "name"]
  plot_title <- paste0("MSC Network ", label_name, " (#", parc, ")")

  plot(
    parc_copy,
    fname = file.path(out_dir, paste0("MSC_", gsub(" ", "_", label_name), ".png")),
    title = plot_title
  )
}
