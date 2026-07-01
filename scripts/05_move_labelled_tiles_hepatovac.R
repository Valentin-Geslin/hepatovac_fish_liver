library(stringr)

working_dir <- "~/Desktop/hepatovac/tiles"
working_dir <- normalizePath(working_dir)  # Convert to absolute path
setwd(working_dir)

move_files_with_pattern <- function(source_folder, destination_folder, pattern) {
  if (!dir.exists(destination_folder)) {
    dir.create(destination_folder, recursive = TRUE)
  }
  
  files_to_move <- list.files(source_folder, pattern = "-labelled", full.names = TRUE)
  
  for (file_path in files_to_move) {
    file_name <- basename(file_path)
    destination_path <- file.path(destination_folder, file_name)
    file.rename(file_path, destination_path)
    message(sprintf("Moved: %s -> %s", file_path, destination_path))
  }
}

subfolders <- list.files(working_dir, full.names = TRUE)

for (subfolder in subfolders) {
  original_dir <- file.path(subfolder, "original_tiles")
  labelled_dir <- file.path(subfolder, "labelled_tiles")
  pattern <- "-labelled"
  
  if (dir.exists(original_dir)) {
    move_files_with_pattern(original_dir, labelled_dir, pattern)
    labelled_files <- list.files(labelled_dir, pattern = pattern, full.names = TRUE)
    new_names <- str_replace(labelled_files, pattern, "")
    file.rename(labelled_files, new_names)
    flsA <- list.files(original_dir, ".png$", full.names = TRUE, recursive = FALSE) 
    flsB <- list.files(labelled_dir, ".png$", full.names = TRUE, recursive = FALSE) 
    ix <- basename(flsB) %in% basename(flsA)
    missing_fls <- flsB[!ix]
    file.remove(missing_fls)
    message(sprintf("Finished processing: %s", subfolder))
  } else {
    message(sprintf("Skipped (missing original_tiles): %s", subfolder))
  }
}
