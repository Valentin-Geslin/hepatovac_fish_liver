library(data.table)
library(tibble)
library(dplyr)
library(readr)

mask_folder <- "~/Desktop/hepatovac/results/masks/"
folders <- list.dirs(mask_folder, recursive = FALSE)
mask_area_list <- vector("list", length(folders))

for (i in seq_along(folders)) {
  folder <- folders[i]
  folder_name <- basename(folder)
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  num_csv_files <- length(csv_files)
  
  if (num_csv_files == 0) {
    mask_area_list[[i]] <- tibble(folder_name, total_area = 0, num_csv_files = 0)
    next
  }
  
  total_area <- 0  
  for (file in csv_files) {
    mask_data <- fread(file, select = "Area", colClasses = list(numeric = "Area")) 
    total_area <- total_area + sum(mask_data$Area, na.rm = TRUE)
  }
  
  mask_area_list[[i]] <- tibble(folder_name, total_area, num_csv_files)
}

mask_hepatovac_area <- bind_rows(mask_area_list)

names(mask_hepatovac_area) <- c("image","hepatovac_area_µm2","num_tiles")
write_csv(mask_hepatovac_area, "~/Desktop/hepatovac/results/mask_hepatovac_area.csv")
