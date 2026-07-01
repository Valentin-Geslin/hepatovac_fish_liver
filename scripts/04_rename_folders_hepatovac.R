working_dir <- "~/Desktop/hepatovac/tiles"
setwd(working_dir)

old_names <- c("TISSUE")
new_names <- c("original_tiles")
subfolders <- list.dirs(working_dir, recursive = FALSE)

for (subfolder in subfolders) {
  for (i in seq_along(old_names)) {
    old_path <- file.path(subfolder, old_names[i])
    new_path <- file.path(subfolder, new_names[i])
    
    if (dir.exists(old_path)) {
      file.rename(old_path, new_path)
      cat("Renamed in", subfolder, ":", old_names[i], "to", new_names[i], "\n")
      extra_folders <- c("labelled_tiles","mask_tiles", "masked_tiles", "selected_tiles")
      for (folder in extra_folders) {
        dir.create(file.path(subfolder, folder), showWarnings = FALSE)
      }
    } else {
      cat("Folder not found in", subfolder, ":", old_names[i], "\n")
    }
  }
}

cat("Folder renaming and creation of additional folders in subfolders completed!\n")
