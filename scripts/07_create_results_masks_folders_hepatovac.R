setwd("~/Desktop/hepatovac/results/masks/")
folders <- list.dirs("~/Desktop/hepatovac/tiles/", recursive = FALSE)

for (folder in folders) {
  folder_name <- basename(folder)
  dir.create(file.path(getwd(), folder_name), showWarnings = FALSE)
}
