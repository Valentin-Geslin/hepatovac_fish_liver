## Load library ------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(bkmr,brms,car,dplyr,GGally,ggplot2,ggpubr,ggcorrplot,ggrepel,glmnet,lme4,mgcv,mice,mixtools,paletteer,readxl,readr,reshape2,stringr,tidyr,tidyverse,purrr)

library(bkmr)
library(brms)
library(car)
library(dplyr)
library(forcats)
library(GGally)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
library(glmnet)
library(lme4)
library(mgcv)
library(mice)
library(paletteer)
library(purrr)
library(readxl)
library(readr)
library(reshape2)
library(stringr)
library(tidyr)
library(purrr)

set.seed(123)

## Import & prepare data -------------------------------------------------------------
setwd("/home/valentin/Desktop/hepatovac/")

### Import image metadata ---------------------------------------------------
image_metadata <- read.table("/home/valentin/Desktop/hepatovac/results/wsi_metadata.txt", row.names=NULL, quote="\"", comment.char="")
image_metadata <- image_metadata[, -c(1,2,3,4,6,8,11,14,17,20,23,25)]
names(image_metadata) <- c("image","","pixel_width","","pixel_height","","resolution","","sizeC","","pixel_width_microns","","pixel_height_microns")
image_metadata <- image_metadata[, -c(2,4,6,8,10,12)]
image_metadata$image <- str_replace_all(image_metadata$image, '.svs', '')
image_metadata$pixel_width <- str_replace_all(image_metadata$pixel_width, ',', '')
image_metadata$pixel_height <- str_replace_all(image_metadata$pixel_height, ',', '')
image_metadata$resolution <- str_replace_all(image_metadata$resolution, ',', '')
image_metadata$sizeC <- str_replace_all(image_metadata$sizeC, ',', '')
image_metadata$pixel_width_microns <- str_replace_all(image_metadata$pixel_width_microns, ',', '')
image_metadata$pixel_width <- as.numeric(image_metadata$pixel_width)
image_metadata$pixel_height <- as.numeric(image_metadata$pixel_height)
image_metadata$resolution <- as.numeric(image_metadata$resolution)
image_metadata$sizeC <- as.numeric(image_metadata$sizeC)
image_metadata$pixel_width_microns <- as.numeric(image_metadata$pixel_width_microns)
image_metadata$image_area_µm2 <- image_metadata$pixel_width*image_metadata$pixel_height*image_metadata$pixel_width_microns*image_metadata$pixel_height_microns

### Import tissue_area ------------------------------------------------------
tissue_area <- read_csv("results/tissue_area_hepatovac.csv")
tissue_area[ , c('Object ID','Object type','Name',"Classification","Parent","ROI","Centroid X µm","Centroid Y µm")] <- list(NULL)
colnames(tissue_area)[colnames(tissue_area) == 'Image'] <- 'image'
colnames(tissue_area)[colnames(tissue_area) == 'Area µm^2'] <- 'total_tissue_area_µm2'
colnames(tissue_area)[colnames(tissue_area) == 'Perimeter µm'] <- 'tissue_perimeter_µm'
tissue_area$image <- str_replace_all(tissue_area$image, '.svs', '')

### Import mask_hepatovac_area --------------------------------------------------------
# Use compile_hepatovac_mask_data.R to compile hepatovac masks csv files 
mask_hepatovac_area <- read.csv("/home/valentin/Desktop/hepatovac/results/mask_hepatovac_area.csv")
names(mask_hepatovac_area) <- c("image","hepatovac_area_µm2","num_tiles")
mask_hepatovac_area$selected_tissue_area_µm2 <- mask_hepatovac_area$num_tiles*422.91^2
mask_hepatovac_area$prop_hepatovac_area <- mask_hepatovac_area$hepatovac_area_µm2*100/mask_hepatovac_area$selected_tissue_area_µm2

### Import biomarker_data ------------------------------------------------------
GIPSA_SEINE_2019 <- read_excel("/home/valentin/Desktop/hepatovac/results/GIPSA_SEINE_2019.xlsx", 
                               col_types = c("text","text", "text", "date", "numeric",
                                             "text", "numeric", "numeric", "text", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "text", "text", "text", 
                                             "text", "text", "text"))
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'individu'] <- 'image'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'taille_cm'] <- 'length_cm'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'poids_g'] <- 'weight_g'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'sexe'] <- 'sex'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'poids_gonades_g'] <- 'gonad_weight_g'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'IGS'] <- 'GSI'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'poids_foie_g'] <- 'liver_weight_g'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'IHP'] <- 'HSI'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'VTG_ng/mL'] <- 'VTG_ng_ml'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'VTG_M/ET'] <- 'VTG_M_ET'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'micronoyaux_/_1000_cellules'] <- 'micronuclei_1000_cells'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'EROD_M/ET'] <- 'EROD_M_ET'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'EROD_pmol/min/mg_prot'] <- 'EROD_pmol_min_mg_prot'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'AChE_M/ET'] <- 'AChE_M_ET'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'AChE_nmol/min/mg_prot'] <- 'AchE_µmol_min_mg_prot'
colnames(GIPSA_SEINE_2019)[colnames(GIPSA_SEINE_2019) == 'COMET_percent_DNA_Tail'] <- 'comet_percentage_DNA_tail'

HQFISH <- read_excel("/home/valentin/Desktop/hepatovac/results/HQFISH.xlsx", 
                          col_types = c("text", "date", "numeric", "text", 
                                        "text", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "text", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "text", "text", 
                                        "text", "text", "text", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric", "numeric", 
                                        "numeric", "numeric"))
colnames(HQFISH)[colnames(HQFISH) == 'code'] <- 'image'

### Merge data frames -----------------------------------------------
data_list <- list(GIPSA_SEINE_2019, HQFISH, image_metadata, mask_hepatovac_area, tissue_area)
data_TOXEM <- data_list %>% reduce(full_join, by='image')
data_TOXEM$image <- str_replace(data_TOXEM$image, "2019", "")
data_TOXEM$image <- str_replace(data_TOXEM$image, "2018", "")

process_data_TOXEM <- function(data_TOXEM) {
  x_columns <- names(data_TOXEM)[grepl("\\.x$", names(data_TOXEM))]
  y_columns <- names(data_TOXEM)[grepl("\\.y$", names(data_TOXEM))]
  for (i in seq_along(x_columns)) {
    new_column_name <- gsub("\\.x$", "", x_columns[i])
    data_TOXEM[[new_column_name]] <- str_replace_all(paste(data_TOXEM[[x_columns[i]]], data_TOXEM[[y_columns[i]]], sep = ""), 'NA', '')
  }
  data_TOXEM <- data_TOXEM %>% dplyr::select(-one_of(c(x_columns, y_columns)))
  
  return(data_TOXEM)
}

data_TOXEM <- process_data_TOXEM(data_TOXEM)

data_TOXEM <- data_TOXEM[, -c(2,3,4)]
#data_TOXEM <- data_TOXEM %>% filter(!is.na(sex))

# Remove columns with NAs
data_TOXEM <- data_TOXEM[,colSums(is.na(data_TOXEM))<nrow(data_TOXEM)]

# Change column type as.numeric
data_TOXEM$hepatovac_area_µm2 <- as.numeric(data_TOXEM$hepatovac_area_µm2)
data_TOXEM$selected_tissue_area_µm2 <- as.numeric(data_TOXEM$selected_tissue_area_µm2)
data_TOXEM$total_tissue_area_µm2 <- as.numeric(data_TOXEM$total_tissue_area_µm2)
data_TOXEM$prop_hepatovac_area <- as.numeric(data_TOXEM$prop_hepatovac_area)

# Rename columns
#data_TOXEM <- data_TOXEM %>% rename(VTG_ng_ml = VTG_ng_mL)
data_TOXEM <- data_TOXEM %>% rename(eleven_KT = `11KT`)
data_TOXEM <- data_TOXEM %>% rename(two_methylnaphtalene_liver_ng_g = `2-methylnaphtalene_liver_ng_g`)
data_TOXEM <- data_TOXEM %>% rename(one_methylnaphtalene_liver_ng_g = `1-methylnaphtalene_liver_ng_g`)
data_TOXEM <- data_TOXEM %>% rename(benzo_e_pyrene_liver_ng_g = `benzo-e-pyrene_liver_ng_g`)
data_TOXEM <- data_TOXEM %>% rename(benzo_ghi_perylene_liver_ng_g = `benzo-ghi-perylene_liver_ng_g`)
data_TOXEM <- data_TOXEM %>% rename(PCB_8 = `PCB-8`)
data_TOXEM <- data_TOXEM %>% rename(PCB_18 = `PCB-18`)
data_TOXEM <- data_TOXEM %>% rename(PCB_28 = `PCB-28`)
data_TOXEM <- data_TOXEM %>% rename(PCB_31 = `PCB-31`)
data_TOXEM <- data_TOXEM %>% rename(PCB_44 = `PCB-44`)
data_TOXEM <- data_TOXEM %>% rename(PCB_49 = `PCB-49`)
data_TOXEM <- data_TOXEM %>% rename(PCB_52 = `PCB-52`)
data_TOXEM <- data_TOXEM %>% rename(PCB_77 = `PCB-77`)
data_TOXEM <- data_TOXEM %>% rename(PCB_101 = `PCB-101`)
data_TOXEM <- data_TOXEM %>% rename(lipids_PI_µg_mg = lipids_PI__µg_mg)
data_TOXEM <- data_TOXEM %>% rename(PCB_110 = `PCB-110`)
data_TOXEM <- data_TOXEM %>% rename(PCB_118 = `PCB-118`)
data_TOXEM <- data_TOXEM %>% rename(PCB_128 = `PCB-128`)
data_TOXEM <- data_TOXEM %>% rename(PCB_132 = `PCB-132`)
data_TOXEM <- data_TOXEM %>% rename(PCB_135 = `PCB-135`)
data_TOXEM <- data_TOXEM %>% rename(PCB_138 = `PCB-138`)
data_TOXEM <- data_TOXEM %>% rename(PCB_149 = `PCB-149`)
data_TOXEM <- data_TOXEM %>% rename(PCB_153 = `PCB-153`)
data_TOXEM <- data_TOXEM %>% rename(PCB_156 = `PCB-156`)
data_TOXEM <- data_TOXEM %>% rename(PCB_169 = `PCB-169`)
data_TOXEM <- data_TOXEM %>% rename(PCB_170 = `PCB-170`)
data_TOXEM <- data_TOXEM %>% rename(PCB_180 = `PCB-180`)
data_TOXEM <- data_TOXEM %>% rename(PCB_187 = `PCB-187`)
data_TOXEM <- data_TOXEM %>% rename(PCB_194 = `PCB-194`)
data_TOXEM <- data_TOXEM %>% rename(endosulfan_alpha = `endosulfan-alpha`)
colnames(data_TOXEM)[99] <- "two_four_dde"
colnames(data_TOXEM)[100] <- "four_four_dde"
colnames(data_TOXEM)[101] <- "two_four_ddd"
colnames(data_TOXEM)[102] <- "four_four_ddd"

rm(list = ls()[! ls() %in% c("data_TOXEM")])

# Recode date as month (categorical variable)
data_TOXEM$date <- as.Date(data_TOXEM$date, format="%Y-%md")
data_TOXEM$month <- str_replace_all(data_TOXEM$, 'January', 'Winter')


data_TOXEM$month <- as.factor(data_TOXEM$month)
data_TOXEM$month <- str_replace_all(data_TOXEM$month, '01', 'January')
data_TOXEM$month <- str_replace_all(data_TOXEM$month, '07', 'July')
data_TOXEM$month <- str_replace_all(data_TOXEM$month, '09', 'September')

# Recode month as season (categorical variable)
data_TOXEM$season <- str_replace_all(data_TOXEM$month, 'January', 'Winter')
data_TOXEM$season <- str_replace_all(data_TOXEM$season, 'July', 'Summer')
data_TOXEM$season <- str_replace_all(data_TOXEM$season, 'September', 'Summer')

# Recode descriptive variables as factors
data_TOXEM$non_specific_lesions_1[data_TOXEM$non_specific_lesions_1 == "NA"] <- NA
data_TOXEM$non_specific_lesions_1[data_TOXEM$non_specific_lesions_1 == "Increased numbers and size of macrophage aggregates"] <- "MMA"
data_TOXEM$prevalence_non_specific_lesions_1 <- ifelse(!is.na(data_TOXEM$non_specific_lesions_1), 1, 0)

data_TOXEM$non_specific_lesions_2[data_TOXEM$non_specific_lesions_2 == "NA"] <- NA
data_TOXEM$prevalence_non_specific_lesions_2 <- ifelse(!is.na(data_TOXEM$non_specific_lesions_2), 1, 0)

data_TOXEM$non_specific_lesions_3[data_TOXEM$non_specific_lesions_3 == "NA"] <- NA
data_TOXEM$non_specific_lesions_3[data_TOXEM$non_specific_lesions_3== "Increased numbers and size of macrophage aggregates"] <- "MMA"
data_TOXEM$prevalence_non_specific_lesions_3 <- ifelse(!is.na(data_TOXEM$non_specific_lesions_3), 1, 0)

data_TOXEM$early_toxicopathic_non_neoplastic_lesions[data_TOXEM$early_toxicopathic_non_neoplastic_lesions == "NA"] <- NA
data_TOXEM$prevalence_early_toxicopathic_non_neoplastic_lesions <- ifelse(!is.na(data_TOXEM$early_toxicopathic_non_neoplastic_lesions), 1, 0)

data_TOXEM$FCA <- paste(data_TOXEM$foci_of_cellular_altreration_FCA, data_TOXEM$FCA, sep="_")
data_TOXEM <- data_TOXEM %>% select(-foci_of_cellular_altreration_FCA)
data_TOXEM$FCA[data_TOXEM$FCA == "NA_NA"] <- NA
data_TOXEM$FCA <- str_replace_all(data_TOXEM$FCA, 'NA_', '')
data_TOXEM$FCA <- str_replace_all(data_TOXEM$FCA, '_NA', '')
data_TOXEM$prevalence_FCA <- ifelse(!is.na(data_TOXEM$FCA), 1, 0)

data_TOXEM <- data_TOXEM %>% rename(TG_FS_ratio = lipids_TG_FS_µg_mg)

# Remove columns not to be include in the analysis
data_TOXEM %>%
   summarise(across(everything(), ~ sum(!is.na(.))))

drop_cols <- c('date','comments', 'resolution', 'sizeC', 'pixel_width', 'pixel_width_microns', 'pixel_height',
               'pixel_height_microns', 'image_area_µm2', 'tissue_perimeter_µm','IHP',
               'sexual_maturation_code_CIEM', 'sexual_maturation_1', 'sexual_maturation_2', 
               'sexual_maturation_3','sexual_maturation_4','sexual_maturation_degeneration', 'sexual_maturation_index',
               'EROD_M_ET', 'AChE_M_ET', 'VTG_M_ET', 'micronuclei_1000_cells')


data_TOXEM <- data_TOXEM %>%
  dplyr::select(-all_of(drop_cols))

data_TOXEM$station <- as.factor(data_TOXEM$station)
data_TOXEM$sex <- as.factor(data_TOXEM$sex)
data_TOXEM$month <- as.factor(data_TOXEM$month)
data_TOXEM$season <- as.factor(data_TOXEM$season)
data_TOXEM$num_tiles <- as.integer(data_TOXEM$num_tiles)
data_TOXEM$gravity_note <- as.factor(data_TOXEM$gravity_note) 
data_TOXEM$sexual_maturation_total <- as.integer(data_TOXEM$sexual_maturation_total)
data_TOXEM$length_cm <- as.numeric(data_TOXEM$length_cm)
data_TOXEM$liver_weight_g <- as.numeric(data_TOXEM$liver_weight_g)
data_TOXEM$weight_g <- as.numeric(data_TOXEM$weight_g)
data_TOXEM$gonad_weight_g <- as.numeric(data_TOXEM$gonad_weight_g)
data_TOXEM$benign_neoplasms <- NA
data_TOXEM$malignant_neoplasms <- NA

# Calculate age_month with Von Bertalanffy growth function
# Y total length (mm) vs X age (month)
L_inf = 424.775
#L_inf = max(data_TOXEM$length_mm)
k = 0.023761
t0 = -2.076
data_TOXEM$length_mm <- data_TOXEM$length_cm*10 
data_TOXEM$age_month <- t0 - (1/k) * log(1 - data_TOXEM$length_mm/L_inf) # Von Bertalanffy growth equation
data_TOXEM$age_year <- data_TOXEM$age_month/12
data_TOXEM$age_year <- round(data_TOXEM$age_year)

data_TOXEM$HSI <- (data_TOXEM$liver_weight_g  * 100) / data_TOXEM$weight_g
data_TOXEM$GSI <- (data_TOXEM$gonad_weight_g * 100) / data_TOXEM$weight_g
data_TOXEM$K_index = (data_TOXEM$length_cm)^3 / data_TOXEM$weight_g
epsilon <- 0.001  
data_TOXEM$log_prop_hepatovac_area <- log(data_TOXEM$prop_hepatovac_area + epsilon)
#data_TOXEM$log_prop_hepatovac_area <- log10(data_TOXEM$prop_hepatovac_area)
data_TOXEM$prop_selected_tissue_area <- data_TOXEM$selected_tissue_area_µm2*100/data_TOXEM$total_tissue_area_µm2

# Reorder columns
data_TOXEM <- data_TOXEM %>% relocate(station, .after = image)
data_TOXEM <- data_TOXEM %>% relocate(sex, .after = station)
data_TOXEM <- data_TOXEM %>% relocate(season, .after = sex)
data_TOXEM <- data_TOXEM %>% relocate(month, .after = season)
data_TOXEM <- data_TOXEM %>% relocate(log_prop_hepatovac_area, .after = month)
data_TOXEM <- data_TOXEM %>% relocate(prop_hepatovac_area, .after = log_prop_hepatovac_area)
data_TOXEM <- data_TOXEM %>% relocate(hepatovac_area_µm2, .after = prop_hepatovac_area)
data_TOXEM <- data_TOXEM %>% relocate(selected_tissue_area_µm2, .after = hepatovac_area_µm2)
data_TOXEM <- data_TOXEM %>% relocate(prop_selected_tissue_area, .after = selected_tissue_area_µm2)
data_TOXEM <- data_TOXEM %>% relocate(total_tissue_area_µm2, .after = prop_selected_tissue_area)
data_TOXEM <- data_TOXEM %>% relocate(num_tiles, .after = total_tissue_area_µm2)

# Fitness index
data_TOXEM <- data_TOXEM %>% relocate(length_cm, .after = num_tiles)
data_TOXEM <- data_TOXEM %>% relocate(length_mm, .after = length_cm)
data_TOXEM <- data_TOXEM %>% relocate(age_year, .after = length_mm)
data_TOXEM <- data_TOXEM %>% relocate(age_month, .after = age_year)
data_TOXEM <- data_TOXEM %>% relocate(K_index, .after = age_month)
data_TOXEM <- data_TOXEM %>% relocate(GR_mm_year, .after = K_index)
data_TOXEM <- data_TOXEM %>% relocate(weight_g, .after = GR_mm_year)
data_TOXEM <- data_TOXEM %>% relocate(carcasse_weight_g, .after = weight_g)
data_TOXEM <- data_TOXEM %>% relocate(liver_weight_g, .after = carcasse_weight_g)
data_TOXEM <- data_TOXEM %>% relocate(HSI, .after = liver_weight_g)

# Reproductive markers
data_TOXEM <- data_TOXEM %>% relocate(gonad_weight_g, .after = HSI)
data_TOXEM <- data_TOXEM %>% relocate(GSI, .after = gonad_weight_g)
data_TOXEM <- data_TOXEM %>% relocate(sexual_maturation_total, .after = GSI)
data_TOXEM <- data_TOXEM %>% relocate(VTG_ng_ml, .after = sexual_maturation_total)
data_TOXEM <- data_TOXEM %>% relocate(E2, .after = VTG_ng_ml)
data_TOXEM <- data_TOXEM %>% relocate(eleven_KT, .after = E2)

# Liver histopathology
data_TOXEM <- data_TOXEM %>% relocate(non_specific_lesions_1, .after = eleven_KT)
data_TOXEM <- data_TOXEM %>% relocate(non_specific_lesions_2, .after = non_specific_lesions_1)
data_TOXEM <- data_TOXEM %>% relocate(non_specific_lesions_3, .after = non_specific_lesions_2)
data_TOXEM <- data_TOXEM %>% relocate(early_toxicopathic_non_neoplastic_lesions, .after = non_specific_lesions_3)
data_TOXEM <- data_TOXEM %>% relocate(FCA, .after = early_toxicopathic_non_neoplastic_lesions)
data_TOXEM <- data_TOXEM %>% relocate(benign_neoplasms, .after = FCA)
data_TOXEM <- data_TOXEM %>% relocate(malignant_neoplasms, .after = benign_neoplasms)
data_TOXEM <- data_TOXEM %>% relocate(gravity_note, .after = malignant_neoplasms)

# Biomarkers
data_TOXEM <- data_TOXEM %>% relocate(TG_FS_ratio, .after = gravity_note)
data_TOXEM <- data_TOXEM %>% relocate(AchE_µmol_min_mg_prot, .after = TG_FS_ratio)
data_TOXEM <- data_TOXEM %>% relocate(EROD_pmol_min_mg_prot, .after = AchE_µmol_min_mg_prot)
data_TOXEM <- data_TOXEM %>% relocate(comet_percentage_DNA_tail, .after = EROD_pmol_min_mg_prot)
data_TOXEM <- data_TOXEM %>% relocate(COMET_OTM, .after = comet_percentage_DNA_tail)
data_TOXEM <- data_TOXEM %>% relocate(TBARS_nmol_eq_MDA_mg_prot, .after = COMET_OTM)

data_TOXEM <- data_TOXEM %>% relocate(prevalence_non_specific_lesions_1, .after = non_specific_lesions_1)
data_TOXEM <- data_TOXEM %>% relocate(prevalence_non_specific_lesions_2, .after = non_specific_lesions_2)
data_TOXEM <- data_TOXEM %>% relocate(prevalence_non_specific_lesions_3, .after = non_specific_lesions_3)
data_TOXEM <- data_TOXEM %>% relocate(prevalence_early_toxicopathic_non_neoplastic_lesions, .after = early_toxicopathic_non_neoplastic_lesions)
data_TOXEM <- data_TOXEM %>% relocate(prevalence_FCA, .after = FCA)

# Recalculate total_lipids
# Neutral lipid (SE, GE, TG, FFA, ALC, FS)
data_TOXEM$total_neutral_lipids_µg_mg <- rowSums(data_TOXEM[ , c("lipids_SE_µg_mg", 
                                                                 "lipids_GE_µg_mg", 
                                                                 "lipids_TG_µg_mg",
                                                                 "lipids_FFA_µg_mg",
                                                                 "lipids_ALC_µg_mg",
                                                                 "lipids_FS_µg_mg")], na.rm=TRUE)

# Storage lipids (SE, GE, TG, FFA, ALC)
data_TOXEM$total_reserve_lipids_µg_mg <- rowSums(data_TOXEM[ , c("lipids_SE_µg_mg", 
                                                                 "lipids_GE_µg_mg", 
                                                                 "lipids_TG_µg_mg",
                                                                 "lipids_FFA_µg_mg",
                                                                 "lipids_ALC_µg_mg")], na.rm=TRUE)

# Polar lipid (SPG, LPC, PC, PS, PI, CL, PE)
data_TOXEM$total_membrane_lipids_µg_mg <- rowSums(data_TOXEM[ , c("lipids_SPG_µg_mg", 
                                                                  "lipids_LPC_µg_mg", 
                                                                  "lipids_PC_µg_mg",
                                                                  "lipids_PS_µg_mg",
                                                                  "lipids_PI_µg_mg",
                                                                  "lipids_CL_µg_mg",
                                                                  "lipids_PE_µg_mg",
                                                                  "lipids_FS_µg_mg")], na.rm=TRUE)

# Membrane lipids (SPG, LPC, PC, PS, PI, CL, PE, FS)
data_TOXEM$total_membrane_lipids_µg_mg <- rowSums(data_TOXEM[ , c("lipids_SPG_µg_mg", 
                                                                 "lipids_LPC_µg_mg", 
                                                                 "lipids_PC_µg_mg",
                                                                 "lipids_PS_µg_mg",
                                                                 "lipids_PI_µg_mg",
                                                                 "lipids_CL_µg_mg",
                                                                 "lipids_PE_µg_mg",
                                                                 "lipids_FS_µg_mg")], na.rm=TRUE)


# Recalculate total PCB
data_TOXEM$total_PCB_liver_mg_kg <- rowSums(select(data_TOXEM, starts_with("PCB_")), na.rm = TRUE)
data_TOXEM$total_PCB_liver_mg_kg[data_TOXEM$total_PCB_liver_mg_kg == 0] <- NA

# Recalculate total_HAP_liver_mg_kg
data_TOXEM$total_HAP_liver_mg_kg <- rowSums(data_TOXEM[ , c("naphtalene_liver_ng_g", 
                                                            "dibenzothiophene_liver_ng_g", 
                                                            "fluoranthene_liver_ng_g",
                                                            "pyrene_liver_ng_g",
                                                            "chrysene_liver_ng_g")], na.rm=TRUE)

# Recalculate total_pesticides_liver_mg_kg
data_TOXEM$total_pesticides_liver_mg_kg <- rowSums(data_TOXEM[ ,c("four_four_dde",
                                                                  "four_four_ddd")], na.rm=TRUE)

data_TOXEM$total_pesticides_liver_mg_kg[data_TOXEM$total_pesticides_liver_mg_kg == 0] <- NA

# Recalculate total_metals_mg_kg
data_TOXEM$total_metals_mg_kg <- rowSums(data_TOXEM[ , c("Ag_liver_mg_kg",
                                                       "Al_liver_mg_kg",
                                                       "Cd_liver_mg_kg",
                                                       "Co_liver_mg_kg",
                                                       "Cr_liver_mg_kg",
                                                       "Cu_liver_mg_kg",
                                                       "Fe_liver_mg_kg",
                                                       "Hg_liver_mg_kg",
                                                       "Mg_liver_mg_kg",
                                                       "Mn_liver_mg_kg",
                                                       "Mo_liver_mg_kg",
                                                       "Ni_liver_mg_kg",
                                                       "Pb_liver_mg_kg",
                                                       "Sr_liver_mg_kg",
                                                       "Ti_liver_mg_kg",
                                                       "V_liver_mg_kg",
                                                       "Zn_liver_mg_kg")], na.rm=TRUE)

data_TOXEM$total_metals_mg_kg[data_TOXEM$total_metals_mg_kg == 0] <- NA

rm(list = ls()[! ls() %in% c("data_TOXEM")])

## Data analysis -----------------------------------------------------------
### prop_hepatovac_area ------------------------------------------------------
summary(data_TOXEM$prop_hepatovac_area)

mean(data_TOXEM$prop_hepatovac_area)
sd(data_TOXEM$prop_hepatovac_area)
median(data_TOXEM$prop_hepatovac_area)
min(data_TOXEM$prop_hepatovac_area)
max(data_TOXEM$prop_hepatovac_area)

data_TOXEM %>%
  group_by(station) %>%
  summarise(mean_prop_hepatovac_area = mean(prop_hepatovac_area), sd_age_month = sd(prop_hepatovac_area))

# prop_hepatovac_area test for normality -> Shapiro test
shapiro.test(data_TOXEM$prop_hepatovac_area) # p-value < 0.05, not normally distributed so we will use log transform variable
wilcox.test(prop_hepatovac_area~station, data = data_TOXEM)

# prop_hepatovac_area ~ sex
ggplot(data_TOXEM) +
  aes(x = sex, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()+
  stat_compare_means(label =  "p.signif", label.x = 1.5)

# prop_hepatovac_area ~ station
ggplot(data_TOXEM) +
  aes(x = station, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()+
  stat_compare_means(label =  "p.signif", label.x = 1.5)

# prop_hepatovac_area ~ month
ggplot(data_TOXEM) +
  aes(x = month, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()+
  stat_compare_means(label =  "p.signif", label.x = 1.5)

# prop_hepatovac_area ~ season
ggplot(data_TOXEM) +
  aes(x = season, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()+
  stat_compare_means(label =  "p.signif", label.x = 1.5)

# prop_hepatovac_area ~ station::season
ggplot(data_TOXEM) +
  aes(x = season, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ station)

ggplot(data_TOXEM) +
  aes(x = season, y = prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  facet_wrap(~ station)+
  stat_compare_means(label =  "p.signif", label.x = 1.5)

# prop_hepatovac_area ~ station::season::sex
ggplot(data_TOXEM) +
  aes(x = season, y = prop_hepatovac_area, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~ station) +
  stat_compare_means(label =  "p.signif", label.x = 1.5) +
  label()

ggplot(data_TOXEM) +
  aes(x = station, y = prop_hepatovac_area, fill = sex) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  facet_wrap(~season) +
  stat_compare_means(aes(group = station), label ="p.signif", label.x = 1.5) +
  labs(x = "Stations", y = "Proportion of tissue affected by hepatocellular vacuolization area (%)")

### log_prop_hepatovac_area -------------------------------------------------

shapiro.test(data_TOXEM$log_prop_hepatovac_area)
qqnorm(data_TOXEM$log_prop_hepatovac_area)
qqline(data_TOXEM$log_prop_hepatovac_area)
hist(data_TOXEM$log_prop_hepatovac_area)
# log_prop_hepatovac_area have a normal distribution so it will be used for statistical analysis

# log_prop_hepatovac_area ~ sex
ggplot(data_TOXEM) +
  aes(x = sex, y = log_prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()

# log_prop_hepatovac_area ~ station
ggplot(data_TOXEM) +
  aes(x = station, y = log_prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()

# log_prop_hepatovac_area ~ month
ggplot(data_TOXEM) +
  aes(x = month, y = log_prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()

# log_prop_hepatovac_area ~ season
ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter()

# log_prop_hepatovac_area ~ station::season
ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  facet_wrap(~ station)

# log_prop_hepatovac_area ~ station::season::sex
ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  facet_wrap(~ station)


### station ~ month ~ sex *** -----------------------------------------------------------
table_station_month_sex <- data_TOXEM %>%
  group_by(station, month, sex) %>%                 # Group by month (season), sex, and station
  summarise(count_samples = n()) %>%                # Count the number of samples
  ungroup()

write_csv(table_station_month_sex, "/home/valentin/Desktop/hepatovac/results/table_station_month_sex.csv")

table(data_TOXEM$month, data_TOXEM$sex, data_TOXEM$station) # unbalanced design (each group have different sample size)
# !!! no sample for Canche in September !!!

leveneTest(log_prop_hepatovac_area ~ sex * station * month, data = data_TOXEM) # Homogeneity of variance among the different groups -> we can proceed with the two-way anova type III

# two-way ANOVA type III (unbalanced design)
aov01 <- aov(log_prop_hepatovac_area ~ sex + station + month, data = data_TOXEM)
Anova(aov01, type = "III")
# Intercept: The intercept represents the mean log-transformed response when all predictors are at their baseline levels. The extremely low p-value indicates that the overall model fit is significant.
# Sex: The F-value of 5.8259 and the p-value of 0.017123 indicate that there is a statistically significant effect of sex on log_prop_hepatovac_area. Since the p-value is less than 0.05, you reject the null hypothesis for sex.
# Station: The F-value of 11.2635 and p-value of 0.001025 suggest a significant effect of station on the response variable. This means different stations have significantly varying impacts on log_prop_hepatovac_area.
# Month: The F-value of 34.6199 and the very small p-value (7.028e-13) indicate a highly significant effect of month on the response. This implies that the month variable can significantly explain the variance in log_prop_hepatovac_area.
# Residuals: The residual sum of squares (4.7894) represents the variation in the response variable not explained by the model. The degree of freedom (136) indicates the amount of data utilized in this analysis.

# Post-hoc test
TukeyHSD(aov01, type = "III")

# sex has a significant effect on log_prop_hepatovac_area
# station does not have a significant effect on log_prop_hepatovac_area. There is no statistically significant difference between Seine and Canche stations for log_prop_hepatovac_area
# month factor significantly affects the response variable, with clear differences between the months examined, particularly highlighting how July and September compare favorably against January.

# Interaction plot
# Interaction plot with faceting by month
ggplot(data_TOXEM, aes(x = station, y = log_prop_hepatovac_area, color = sex)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  facet_wrap(~ month) +
  labs(title = "Interaction Plot: Sex and Station by month",
       x = "Station",
       y = "Mean Log Proportion Hepatovac Area") +
  theme_minimal()

ggplot(data_TOXEM, aes(x = station, y = log_prop_hepatovac_area, color = month)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean", aes(group = month)) +
  labs(title = "Interaction Plot: Sex and Station Colored by month",
       x = "Station",
       y = "Mean Log Proportion Hepatovac Area") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Check ANOVA assumption
hist(aov01$residuals)
qqnorm(aov01$residuals)
qqline(aov01$residuals)
bartlett.test(log_prop_hepatovac_area ~ sex, data=data_TOXEM) # p-value > 0.05 -> equal variance
bartlett.test(log_prop_hepatovac_area ~ station, data=data_TOXEM) #  p-value > 0.05 -> equal variance
bartlett.test(log_prop_hepatovac_area ~ month, data=data_TOXEM) #  p-value > 0.05 -> equal variance

boxplot(log_prop_hepatovac_area ~ sex, data = data_TOXEM)
wilcox.test(log_prop_hepatovac_area ~ sex, data = data_TOXEM) # p-value > 0.05 -> no significant difference between sex
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$sex, p.adjust.method = "bonferroni")

boxplot(log_prop_hepatovac_area ~ station, data = data_TOXEM)
wilcox.test(log_prop_hepatovac_area ~ station, data = data_TOXEM) # p-value > 0.05 ->no significant difference between station
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$station, p.adjust.method = "bonferroni")

boxplot(log_prop_hepatovac_area ~ month, data = data_TOXEM)
kruskal.test(log_prop_hepatovac_area ~ month, data = data_TOXEM) # p-value < 0.05 -> significant difference between month
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$month, p.adjust.method = "bonferroni")

### station ~ season ~ sex *** -----------------------------------------------------------
table_station_season_sex <- data_TOXEM %>%
  group_by(station, season, sex) %>%                 
  summarise(count_samples = n()) %>%
  ungroup()

write_csv(table_station_season_sex, "/home/valentin/Desktop/hepatovac/results/table_station_season_sex.csv")


table(data_TOXEM$season, data_TOXEM$sex, data_TOXEM$station) # unbalanced design (each group have different sample size) -> two-way anova type III (unbalanced design)

leveneTest(log_prop_hepatovac_area ~ station * season * sex, data = data_TOXEM) # p-value > 0.05 -> equal variance among the different groups station * season * sex
leveneTest(log_prop_hepatovac_area ~ station * season, data = data_TOXEM) # p-value < 0.05 -> variance not equal among the different groups station * season 
leveneTest(log_prop_hepatovac_area ~ station, data = data_TOXEM) # p-value > 0.05 -> equal variance among the different stations
# Homogeneity of variance among the different groups -> ANOVA’s are fairly robust to unequal sample sizes if the variances across each treatment combination are still equal

data_TOXEM %>%
  group_by(station, season) %>%
  summarise(mean(prop_hepatovac_area))

# two-way ANOVA type III (unbalanced design)
aov01 <- aov(log_prop_hepatovac_area ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
Anova(aov01, type = "III")

# Intercept: The intercept represents the mean log-transformed response when all predictors are at their baseline levels. The extremely low p-value indicates that the overall model fit is significant.
# Station: The F-value of 9.78 and p-value < 0.05 indicate a significant effect of station on the response variable.
# Season: The F-value of 25.8 and the p-value < 0.05 indicate a significant effect of season on the response. This implies that the season variable can significantly explain the variance in log_prop_hepatovac_area.
# Sex: The F-value of 2.17 and the p-value > 0.05 indicate that there is no statistically significant effect of sex on log_prop_hepatovac_area.
# Residuals: The residual sum of squares (5.099) represents the variation in the response variable not explained by the model. The degree of freedom (133) indicates the amount of data used in this analysis.

# Wilcoxon test for non-parametric data
wilcox.test(log_prop_hepatovac_area ~ station, data = data_TOXEM) # p-value < 0.05 -> significant difference between station
wilcox.test(log_prop_hepatovac_area ~ season, data = data_TOXEM) # p-value < 0.05 -> significant difference between season
wilcox.test(log_prop_hepatovac_area ~ sex, data = data_TOXEM) # p-value > 0.05 -> no significant difference between sex

# Post-hoc test

## Tukey HSD test -> account for unbalanced design
TukeyHSD(aov01, type = "III")
# Statistically significant differences found between mean log_prop_hepatovac_area for the following pairs:

# Seine-Canche 0.08816944 -> significant difference between the two stations on the mean log_prop_hepatovac_area. There is a statistically significant difference between Seine and Canche stations, with individuals from Seine presenting a higher rate of log_prop_hepatovac_area than individuals from the Canche station.

# Winter-Summer -0.1575005 -> significant difference between the two season on the mean log_prop_hepatovac_area. Individuals sampled during winter showing a significantly lower rate of of log_prop_hepatovac_area than individuals sampled in summer.

# sex -> no significant difference between females and males on the rate of log_prop_hepatovac_area.

# Canche:Winter-Canche:Summer -0.15549011 -> significant difference between winter-summer at Canche station, with fish sampled in Summer presenting a higher rate of log_prop_hepatovac_area than fish sampled at Canche in winter
# Canche:Winter-Seine:Summer  -0.23469115 -> significant difference between the two stations between winter and summer, with Canche in winter showing lower rate of log_prop_hepatovac_area than Seine in summer.
# Seine:Winter-Seine:Summer   -0.18671289 -> significant difference between the two seasons at the Seine station, with winter showing lower rate of log_prop_hepatovac_area than summer.

# Seine:M-Canche:F   0.155461881
# Seine:M-Canche:M   0.163461960

# Winter:F-Summer:F -0.18725209 -> significant difference between the two seasons for females, with females in winter showing lower rate of log_prop_hepatovac_area than females in summer.
# Summer:M-Winter:F  0.21371167 -> significant difference between males sampled in summer and females sampled in winter, with the males in summer presenting a higher rate of log_hepatovac_area than females in winter.

# Canche:Winter:F-Canche:Summer:F -0.23649614 -> significant difference between the two seasons for females at the Canche station, with females sampled in winter showing lower rate of log_prop_hepatovac_area than females sampled in summer.
# Canche:Winter:F-Seine:Summer:F  -0.331667533 -> significant difference between the two stations between winter and summer for females, with females in Canche in winter showing lower rate of log_prop_hepatovac_area than females in Seine in summer.
# Seine:Winter:F-Seine:Summer:F   -0.24187317 -> significant difference between the two seasons for females at the Seine station, with females in winter showing lower rate of log_prop_hepatovac_area than females in summer.
# Seine:Summer:M-Canche:Winter:F   0.26414979 -> significant difference between males in Seine in summer and females in Canche in winter, with males in Seine in summer showing higher rate of log_prop_hepatovac_area than females in Canche in winter.

plot(TukeyHSD(aov01, conf.level=.95), las = 2)

ggplot(data_TOXEM) +
  aes(x = station, y = log_prop_hepatovac_area) +
  geom_boxplot() +
  stat_compare_means(label =  "p.signif", label.x = 1.5)

ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area) +
  geom_boxplot() +
  stat_compare_means(label =  "p.signif", label.x = 1.5)

ggplot(data_TOXEM) +
  aes(x = sex, y = log_prop_hepatovac_area) +
  geom_boxplot() +
  stat_compare_means(label =  "p.signif", label.x = 1.5)

ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area) +
  geom_boxplot() +
  facet_wrap(~ station) + stat_compare_means(aes(group = season), label =  "p.signif", label.x = 1.5)

ggplot(data_TOXEM) +
  aes(x = season, y = log_prop_hepatovac_area, fill = sex) +
  geom_boxplot() +
  facet_wrap(~ station) +
  labs(x = "", y = "Mean Log Proportion Hepatovac Area")

## Sheffé test -> more conservative than TukeyHSD (https://www.statology.org/scheffe-test-in-r/)
library(DescTools)
ScheffeTest(aov01)
# Statistically significant differences found between mean log_prop_hepatovac_area for the following pairs:

# Winter-Summer -0.1575005 -> significant difference between the two seasons, with winter showing lower rate of log_prop_hepatovac_area than summer.

# Seine:Winter-Seine:Summer   -0.200804298 -> significant difference between the two seasons at the Seine station, with winter showing lower rate of log_prop_hepatovac_area than summer.
# Canche:Winter-Seine:Summer  -0.204387753 -> significant difference between the two stations between winter and summer, with Canche in winter showing lower rate of log_prop_hepatovac_area than Seine in summer.

# Winter:F-Summer:F -0.22500745 -> significant difference between the two seasons for females, with females in winter showing lower rate of log_prop_hepatovac_area than females in summer

# Seine:Winter:F-Seine:Summer:F   -0.24187317 -> significant difference between the two seasons for females at the Seine station, with females in winter showing lower rate of log_prop_hepatovac_area than females in summer.
# Canche:Winter:F-Seine:Summer:F  -0.331667533 -> significant difference between the two stations between winter and summer for females, with females in Canche in winter showing lower rate of log_prop_hepatovac_area than females in Seine in summer.
# Seine:Summer:M-Canche:Winter:F   0.26414979 -> significant difference between males in Seine in summer and females in Canche in winter, with males in Seine in summer showing higher rate of log_prop_hepatovac_area than females in Canche in winter.


# Interaction plot
ggplot(data_TOXEM, aes(x = season, y = log_prop_hepatovac_area, color = sex)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean") +
  facet_wrap(~ station) +
  labs(title = "Interaction Plot: Station ~ Season ~ Sex",
       x = "Season",
       y = "Mean Log Proportion Hepatovac Area") +
  theme_minimal()

ggplot(data_TOXEM, aes(x = station, y = log_prop_hepatovac_area, color = season)) +
  stat_summary(geom = "point", fun = "mean", size = 3) +
  stat_summary(geom = "line", fun = "mean", aes(group = season)) +
  labs(title = "Interaction Plot: Station ~ Season",
       x = "Station",
       y = "Mean Log Proportion Hepatovac Area") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Check ANOVA assumption
hist(aov01$residuals)
qqnorm(aov01$residuals)
qqline(aov01$residuals)
bartlett.test(log_prop_hepatovac_area ~ station, data=data_TOXEM) #  p-value > 0.05 -> equal variance
bartlett.test(log_prop_hepatovac_area ~ season, data=data_TOXEM) #  p-value > 0.05 -> equal variance
bartlett.test(log_prop_hepatovac_area ~ sex, data=data_TOXEM) # p-value > 0.05 -> equal variance

boxplot(log_prop_hepatovac_area ~ station, data = data_TOXEM)
wilcox.test(log_prop_hepatovac_area ~ station, data = data_TOXEM) # p-value < 0.05 -> significant difference between station
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$station, p.adjust.method = "bonferroni")

boxplot(log_prop_hepatovac_area ~ season, data = data_TOXEM)
wilcox.test(log_prop_hepatovac_area ~ season, data = data_TOXEM) # p-value < 0.05 -> significant difference between season
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$season, p.adjust.method = "bonferroni")

boxplot(log_prop_hepatovac_area ~ sex, data = data_TOXEM)
wilcox.test(log_prop_hepatovac_area ~ sex, data = data_TOXEM) # p-value > 0.05 -> no significant difference between sex
pairwise.wilcox.test(data_TOXEM$log_prop_hepatovac_area, data_TOXEM$sex, p.adjust.method = "bonferroni")

rm(list = ls()[! ls() %in% c("data_TOXEM")])
## total_tissue_area_µm2 ---------------------------------
ggplot(data_TOXEM, aes(x = image, y = total_tissue_area_µm2)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

# log_prop_hepatovac_area ~ total_tissue_area_µm2 
ggplot(data_TOXEM)+ 
  aes(total_tissue_area_µm2, log_prop_hepatovac_area, col = station) +
  geom_point() 
#  geom_smooth(method = lm,  se = F)

cor.test(x=data_TOXEM$total_tissue_area_µm2, y=data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_total_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ total_tissue_area_µm2, data = data_TOXEM)
summary(lm_total_tissue_area_µm2)

poly2_total_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ poly(total_tissue_area_µm2, 2), data = data_TOXEM)
summary(poly2_total_tissue_area_µm2)

poly3_total_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ poly(total_tissue_area_µm2, 3), data = data_TOXEM)
summary(poly3_total_tissue_area_µm2)

gam_total_tissue_area_µm2 <- gam(log_prop_hepatovac_area ~ s(total_tissue_area_µm2), data = data_TOXEM)
summary(gam_total_tissue_area_µm2)

# no statistically significant relationship between response variable and total_tissue_area_µm2

## selected_tissue_area_µm2 ---------------------------------
ggplot(data_TOXEM)+ 
  aes(selected_tissue_area_µm2, log_prop_hepatovac_area) +
  geom_point()

cor.test(x=data_TOXEM$selected_tissue_area_µm2, y=data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_selected_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ selected_tissue_area_µm2, data = data_TOXEM)
summary(lm_selected_tissue_area_µm2)

poly2_selected_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ poly(selected_tissue_area_µm2, 2), data = data_TOXEM)
summary(poly2_selected_tissue_area_µm2)

poly3_selected_tissue_area_µm2 <- lm(log_prop_hepatovac_area ~ poly(selected_tissue_area_µm2, 3), data = data_TOXEM)
summary(poly3_selected_tissue_area_µm2)

gam_selected_tissue_area_µm2 <- gam(log_prop_hepatovac_area ~ s(selected_tissue_area_µm2), data = data_TOXEM)
summary(gam_selected_tissue_area_µm2)

# no statistically significant relationship between response variable and selected_tissue_area_µm2

## prop_selected_tissue_area *** ---------------------------------
mean(data_TOXEM$prop_selected_tissue_area)
median(data_TOXEM$prop_selected_tissue_area)
sd(data_TOXEM$prop_selected_tissue_area)
min(data_TOXEM$prop_selected_tissue_area)
max(data_TOXEM$prop_selected_tissue_area)

ggplot(data_TOXEM, aes(x = image, y = prop_selected_tissue_area)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

# length_cm ~ prop_selected_tissue_area
ggplot(data_TOXEM)+ 
  aes(length_cm, prop_selected_tissue_area) +
  geom_point() +
  geom_smooth(method = lm,  se = T)

cor.test(y=data_TOXEM$prop_selected_tissue_area, x=data_TOXEM$length_cm, method = "pearson", use = "complete.obs")
# statistically significant positive correlation ->  prop_selected_tissue_area tend to increase when length_cm increase

# age_month ~ prop_selected_tissue_area
ggplot(data_TOXEM)+ 
  aes(age_month, prop_selected_tissue_area) +
  geom_point() +
  geom_smooth(method = lm,  se = T)

cor.test(y=data_TOXEM$prop_selected_tissue_area, x=data_TOXEM$age_month, method = "pearson", use = "complete.obs")
# statistically significant positive correlation ->  prop_selected_tissue_area tend to increase when age_month increase

# log_prop_hepatovac_area ~ prop_selected_tissue_area
ggplot(data_TOXEM)+ 
  aes(prop_selected_tissue_area, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = T)

ggplot(data_TOXEM)+ 
  aes(prop_selected_tissue_area, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = NULL,  se = T)

ggplot(data_TOXEM) + 
  aes(prop_selected_tissue_area, log_prop_hepatovac_area, group = season, colour = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F)
#stat_smooth(method='lm', formula = y ~ poly(x,3), size = 1, se = F)

cor.test(x=data_TOXEM$prop_selected_tissue_area, y=data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# statistically significant negative correlation -> log_prop_hepatovac_area tend to decrease when prop_selected_tissue_area increase
# need to integrate prop_selected_tissue into the next models to account for this relationship

# lm_prop_selected_tissue_area <- lm(log_prop_hepatovac_area ~ prop_selected_tissue_area + season, data = data_TOXEM)
# summary(lm_prop_selected_tissue_area)
# 
# poly2_prop_selected_tissue_area <- lm(log_prop_hepatovac_area ~ poly(prop_selected_tissue_area, 2) + season, data = data_TOXEM)
# summary(poly2_prop_selected_tissue_area)
# 
# poly3_prop_selected_tissue_area <- lm(log_prop_hepatovac_area ~ poly(prop_selected_tissue_area, 3) + season, data = data_TOXEM)
# summary(poly3_prop_selected_tissue_area)
# 
# gam_prop_selected_tissue_area <- gam(log_prop_hepatovac_area ~ s(prop_selected_tissue_area) + season, data = data_TOXEM)
# summary(gam_prop_selected_tissue_area)

# AIC(lm_prop_selected_tissue_area, poly2_prop_selected_tissue_area, poly3_prop_selected_tissue_area, gam_prop_selected_tissue_area)

library(data.table)
tissue_area <- data_TOXEM %>% select(image, total_tissue_area_µm2, selected_tissue_area_µm2)
tissue_area <- reshape2::melt(tissue_area, id.vars = "image", variable.name = "tissue", value.name="area_µm2")

ggplot(tissue_area, aes(x = image , y= area_µm2, fill = tissue)) +
  geom_bar(position="dodge", stat = "identity")

rm(tissue_area)

## num_tiles ---------------------------------
ggplot(data_TOXEM, aes(x = image, y = num_tiles)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

# num_tiles ~ log_prop_hepatovac_area
ggplot(data_TOXEM)+ 
  aes(num_tiles, log_prop_hepatovac_area) +
  geom_point() + 
  geom_smooth(method = NULL, se = F)

cor.test(x=data_TOXEM$num_tiles, y=data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_num_tiles <- lm(log_prop_hepatovac_area ~ num_tiles, data = data_TOXEM)
summary(lm_num_tiles)

poly2_num_tiles <- lm(log_prop_hepatovac_area ~ poly(num_tiles, 2), data = data_TOXEM)
summary(poly2_num_tiles)

poly3_num_tiles <- lm(log_prop_hepatovac_area ~ poly(num_tiles, 3), data = data_TOXEM)
summary(poly3_num_tiles)

gam_num_tiles <- gam(log_prop_hepatovac_area ~ s(num_tiles), data = data_TOXEM)
summary(gam_num_tiles)
# no statistically significant relationship between response variable and num_tiles

## fitness_data_TOXEM ----------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])
### length_cm *** --------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = length_cm, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station) +
  stat_compare_means(aes(group = season), label ="p.signif", label.x = 1.5) +
  labs(x = "Stations", y = "Total body length (cm)")

data_TOXEM %>%
  group_by(station, season) %>%
  summarise(mean_length_cm = mean(length_cm), sd_length_cm = sd(length_cm), min_length_cm = min(length_cm), max_length_cm = max(length_cm))

aov01 <- aov(length_cm ~ station + season + sex,  data = data_TOXEM)
Anova(aov01, type = "III")
# Intercept: The intercept represents the mean log-transformed response when all predictors are at their baseline levels. The extremely low p-value indicates that the overall model fit is significant.
# Station: The F-value of 69.89 and p-value < 0.05 suggest a significant effect of station on the response variable. This means different stations have significantly varying impacts on log_prop_hepatovac_area.
# Station: The F-value of 21.47 and p-value < 0.05 suggest a significant effect of season on the response variable. This means different seasons have significantly varying impacts on log_prop_hepatovac_area.
# Sex: The F-value of 3.24 and the p-value > 0.05 indicate no significant effect of sex on log_prop_hepatovac_area. 
# Residuals: The residual sum of squares represents the variation in the response variable not explained by the model. The degree of freedom (137) indicates the amount of data utilized in this analysis.

TukeyHSD(aov01, type = "III")
# Significant difference in length_cm between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_TOXEM)+ 
  aes(length_cm, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_TOXEM$length_cm, y=data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_length_cm <- lm(log_prop_hepatovac_area ~ length_cm + station * season * sex, data = data_TOXEM)
summary(lm_length_cm)
# length_cm does not significantly affect the response in this model

# Non linear regression
# GAM
ggplot(data_TOXEM)+ 
  aes(length_cm, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_length_cm <- gam(log_prop_hepatovac_area ~ s(length_cm) + station * season * sex, data = data_TOXEM)
summary(gam_length_cm)
# length_cm does not significantly affect the response in this model.

# Polynomial
ggplot(data_TOXEM) + 
  aes(length_cm, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_length_cm <- lm(log_prop_hepatovac_area ~ poly(length_cm, 2) + station * season, data = data_TOXEM)
summary(poly2_length_cm)

ggplot(data_TOXEM) + 
  aes(length_cm, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  stat_smooth(method='lm', formula = y ~ poly(x,3), size = 1, se = F) +
#  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ station)

poly3_length_cm <- lm(log_prop_hepatovac_area ~ poly(length_cm, 3) + station * season, data = data_TOXEM)
summary(poly3_length_cm)
# length_cm significantly affect the response in this model.

poly3_length_cm <- lm(log_prop_hepatovac_area ~ poly(length_cm, 3) + season, data = data_TOXEM)
summary(poly3_length_cm)

ggplot(data_TOXEM) +
  aes(length_cm, log_prop_hepatovac_area, col = season, group = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F)

ggplot(data_TOXEM) +
  aes(length_cm, log_prop_hepatovac_area, col = station, group = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ season)

# ggplot(data_TOXEM) +
#   aes(length_cm, log_prop_hepatovac_area, col = season, group = season) +
#   geom_point() +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
#   facet_wrap(~ station)

# Compare models
# vif -> colinearity between variables
# vif < 10 -> no colinearity between variables
vif(poly3_length_cm) # no colinearity between variables

summary(poly3_length_cm)$r.squared

AIC(poly2_length_cm, poly3_length_cm)
# Lower values indicate better model fit, penalizing complexity.

anova(poly2_length_cm, poly3_length_cm, test = "F")
# p-value < 0.05 which means that the gam model provides a significantly better fit to the data than the other model.

# Plot fitted values
plot(data_TOXEM$length_cm, data_TOXEM$log_prop_hepatovac_area, main="Actual vs Fitted")
points(data_TOXEM$length_cm, fitted(lm_length_cm), col="green", pch=4)
points(data_TOXEM$length_cm, fitted(gam_length_cm), col="blue", pch=3)
points(data_TOXEM$length_cm, fitted(poly2_length_cm), col="red", pch=2)
points(data_TOXEM$length_cm, fitted(poly3_length_cm), col="black", pch=5)
legend("topleft", legend=c("LM","GAM","PM"), col=c("green","blue","red"), pch=c(2, 3))

plot(data_TOXEM$length_cm, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$length_cm, predict(lm_length_cm), col = "green")
points(data_TOXEM$length_cm, predict(gam_length_cm), col = "blue")
points(data_TOXEM$length_cm, predict(poly2_length_cm), col = "red")
points(data_TOXEM$length_cm, predict(poly3_length_cm), col="black")

# Check assumptions
plot(poly3_length_cm)
res <- resid(poly3_length_cm)
plot(fitted(poly3_length_cm), res) # homogeneity of the variances
abline(0,0)
qqnorm(res)
qqline(res) 
plot(density(res)) # normality of the variances


### age_month ***----------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(age_month))

ggplot(data_TOXEM) +
  aes(x = season, y = age_month, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_TOXEM %>%
  group_by(season, station, sex) %>%
  summarise(mean_age_month = mean(age_month), sd_age_month = sd(age_month))

aov01 <- aov(age_month ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

(age_month_table <- data_TOXEM %>%
    group_by(station, season, sex) %>%
    summarise_at(vars(age_month), list(mean_age_month = mean, median_age_month = median, min_age_month = min, max_age_month = max)))

write_csv(age_month_table, "/home/valentin/Desktop/hepatovac/results/age_month_table.csv")

library(plyr)

mu <- ddply(data_TOXEM, "station", summarise, grp.mean=mean(age_month))
ggplot(data_TOXEM, aes(x=age_month, color=station)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=station),
             linetype="dashed")

mu <- ddply(data_TOXEM, "season", summarise, grp.mean=mean(age_month))
ggplot(data_TOXEM, aes(x=age_month, color=season)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=season),
             linetype="dashed")

mu <- ddply(data_TOXEM, "sex", summarise, grp.mean=mean(age_month))
ggplot(data_TOXEM, aes(x=age_month, color=sex)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")

ggplot(data_TOXEM) +
  aes(age_month, length_mm, shape = sex) +
  geom_point() +
  geom_text(label=data_TOXEM$image)

# Boxplot per station & season & sex
ggplot(data_TOXEM) +
  aes(season, age_month, fill = sex) +
  geom_boxplot() +
  facet_wrap(~station)

aov01 <- aov(age_month ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in age_month between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_TOXEM$age_month, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_age_month <- lm(log_prop_hepatovac_area ~ age_month + station + season + sex, data = data_TOXEM)
summary(lm_age_month)
# age_month does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_age_month <- gam(log_prop_hepatovac_area ~ s(age_month) + station * season, data = data_TOXEM)
summary(gam_age_month)
# age_month does not significantly affect the response in this model.

# Polynomial
(age_month <- ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station))

poly2_age_month <- lm(log_prop_hepatovac_area ~ poly(age_month, 2) + station * season * sex, data = data_clean)
summary(poly2_age_month)
# age_month significantly affect the response in this model.

ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  stat_smooth(method='lm', formula = y ~ poly(x,3) + season, size = 1, se = F) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ station)

poly3_age_month <- lm(log_prop_hepatovac_area ~ poly(age_month, 3) + station * season * sex, data = data_clean)
summary(poly3_age_month)
# age_month has a non-linear (cubic) significant effect on log_prop_hepatovac_area in this model.

poly3_age_month <- lm(log_prop_hepatovac_area ~ poly(age_month, 3) + season, data = data_TOXEM)
summary(poly3_age_month)

ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F)

ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = station, group = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ season)

ggplot(data_TOXEM) +
  aes(age_month, log_prop_hepatovac_area, col = season, group = season) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ station)

# Compare models
vif(lm_age_month)
vif(gam_age_month)
vif(poly2_age_month)
vif(poly3_age_month)
# vif < 10 -> no collinearity between variables

summary(lm_age_month)$r.squared
summary(gam_age_month)$r.sq
summary(poly2_age_month)$r.squared
summary(poly3_age_month)$r.squared

AIC(poly2_age_month, poly3_age_month)
BIC(poly2_age_month, poly3_age_month)
# Lower values indicate better model fit, penalizing complexity.

anova(poly2_age_month, poly3_age_month, test = "F")

# Plot fitted values
plot(data_TOXEM$age_month, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$age_month, fitted(lm_age_month), col = "green", pch = 4)
points(data_TOXEM$age_month, fitted(gam_age_month), col = "blue", pch = 3)
points(data_TOXEM$age_month, fitted(poly2_age_month), col = "red", pch = 2)
points(data_TOXEM$age_month, fitted(poly3_age_month), col = "orange", pch = 5)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$age_month, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$age_month, predict(lm_age_month), col = "green")
points(data_TOXEM$age_month, predict(gam_age_month), col = "blue")
points(data_TOXEM$age_month, predict(poly2_age_month), col = "red")

# Check assumptions
plot(poly2_age_month)
res <- resid(poly2_age_month)
plot(fitted(poly2_age_month), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly2_age_month model show a significant effect of age_month on log_prop_hepatovac_area.
# lm(log_prop_hepatovac_area ~ variable + station * season * sex + poly(age_month, 2), data = data_TOXEM)

### age_year ---------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = age_year, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~station)

age_year_table <- data_TOXEM %>%
  group_by(station, season, sex) %>%
  summarise_at(
    vars(age_year),
    list(mean_age_year = mean, min_age_year = min, max_age_year = max))

write_csv(age_year_table, "/home/valentin/Desktop/hepatovac/results/age_year_table.csv")

library(plyr)

mu <- ddply(data_TOXEM, "station", summarise, grp.mean=mean(age_year))
ggplot(data_TOXEM, aes(x=age_year, color=station)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=station),
             linetype="dashed")

mu <- ddply(data_TOXEM, "season", summarise, grp.mean=mean(age_year))
ggplot(data_TOXEM, aes(x=age_year, color=season)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=season),
             linetype="dashed")

mu <- ddply(data_TOXEM, "sex", summarise, grp.mean=mean(age_year))
ggplot(data_TOXEM, aes(x=age_year, color=sex)) +
  geom_histogram(fill="white", position="dodge")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")

### K_index ***--------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = K_index, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(K_index ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in K_index between males & females and between fishes sampled in Winter vs Summer.

# Linear regression
(K_index <- ggplot(data_TOXEM) +
  aes(K_index, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station))

cor.test(x = data_TOXEM$K_index, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_K_index <- lm(log_prop_hepatovac_area ~ K_index, data = data_TOXEM)
summary(lm_K_index)
# lm_K_index have a significant effect on log_prop_hepatovac_area in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(K_index, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_K_index <- gam(log_prop_hepatovac_area ~ s(K_index) + station * season * sex, data = data_TOXEM)
summary(gam_K_index)
# gam_K_index shows a significant effect on log_prop_hepatovac_area in this model.

# Polynomial
ggplot() +
  aes(K_index, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_K_index <- lm(log_prop_hepatovac_area ~ poly(K_index, 2) + station * season * sex, data = data_TOXEM)
summary(poly2_K_index)
# poly2_K_index model shows a significant effect on log_prop_hepatovac_area in this model.

ggplot() +
  aes(K_index, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_K_index <- lm(log_prop_hepatovac_area ~ poly(K_index, 3) + station * season * sex, data = data_TOXEM)
summary(poly3_K_index)
# poly3_K_index model shows a significant effect on log_prop_hepatovac_area in this model.

# Compare models
vif(lm_K_index)
vif(gam_K_index)
vif(poly2_K_index)
vif(poly3_K_index)
# vif < 10 -> no collinearity between variables

summary(lm_K_index)$r.squared
summary(gam_K_index)$r.sq
summary(poly2_K_index)$r.squared
summary(poly3_K_index)$r.squared

AIC(lm_K_index, gam_K_index, poly2_K_index, poly3_K_index)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_K_index, gam_K_index, poly2_K_index, poly3_K_index, test = "F")

# Plot fitted values
plot(data_TOXEM$K_index, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points($K_index, fitted(lm_K_index), col = "green", pch = 4)
points($K_index, fitted(gam_K_index), col = "blue", pch = 3)
points($K_index, fitted(poly2_K_index), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$K_index, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$K_index, predict(lm_K_index), col = "green")
points(data_TOXEM$K_index, predict(gam_K_index), col = "blue")
points(data_TOXEM$K_index, predict(poly2_K_index), col = "red")

# Check assumptions
plot(gam_K_index)
res <- resid(gam_K_index)
plot(fitted(gam_K_index), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# gam_K_index model show a significant effect of K_index on log_prop_hepatovac_area.
# gam_K_index

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
 <- na.omit(data_TOXEM[, c("K_index", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot()+ 
  aes(K_index, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=$K_index, y=$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_K_index <- lm(EROD_pmol_min_mg_prot ~ K_index + station * season * sex, data = )
summary(lm_K_index)
# K_index significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_K_index <- gam(EROD_pmol_min_mg_prot ~ s(K_index) + station * season * sex, data = )
summary(gam_K_index)
# K_index significantly affect the response in this model.

# Polynomial
poly2_K_index <- lm(EROD_pmol_min_mg_prot ~ poly(K_index, 2) + station * season * sex, data = )
summary(poly2_K_index)
# K_index significantly affect the response in this model.

poly3_K_index <- lm(EROD_pmol_min_mg_prot ~ poly(K_index, 3) + station * season * sex, data = )
summary(poly3_K_index)
# K_index does not significantly affect the response in this model.

### GR_mm_year --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(GR_mm_year))

ggplot(data_TOXEM) +
  aes(season, GR_mm_year, fill = sex) +
  geom_boxplot()+
  #geom_col(stat = "identity", position = "dodge") +
  facet_wrap(~station)

aov01 <- aov(GR_mm_year ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in GR_mm_year between fishes from Canche-Summer and Seine-Summer -> GR significantly higher in fishes from Seine samples in summer comparatively to fishes sampled in Canche at the same season

data_clean %>%
  group_by(station) %>%
  summarise(mean_GR_mm_year = mean(GR_mm_year),
            sd_GR_mm_year = sd(GR_mm_year))

# Linear regression
ggplot(data_TOXEM) +
  aes(GR_mm_year, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_TOXEM$GR_mm_year, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# not significant

lm_GR_mm_year <- lm(log_prop_hepatovac_area ~ GR_mm_year + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_GR_mm_year)

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(GR_mm_year, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_GR_mm_year <- gam(log_prop_hepatovac_area ~ s(GR_mm_year) + station * season * sex, data = data_TOXEM)
summary(gam_GR_mm_year)

# Polynomial
ggplot(data_TOXEM) +
  aes(GR_mm_year, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_GR_mm_year <- lm(log_prop_hepatovac_area ~ poly(GR_mm_year, 2) + station * season, data = data_TOXEM)
summary(poly2_GR_mm_year)

ggplot(data_TOXEM) +
  aes(GR_mm_year, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_GR_mm_year <- lm(log_prop_hepatovac_area ~ poly(GR_mm_year, 3) + station * season, data = data_TOXEM)
summary(poly3_GR_mm_year)

# Compare models
vif(lm_GR_mm_year)
vif(gam_GR_mm_year)
vif(poly2_GR_mm_year)
vif(poly3_GR_mm_year)
# vif < 10 -> no collinearity between variables

summary(lm_GR_mm_year)$r.squared
summary(gam_GR_mm_year)$r.sq
summary(poly2_GR_mm_year)$r.squared
summary(poly3_GR_mm_year)$r.squared

AIC(lm_GR_mm_year, gam_GR_mm_year, poly2_GR_mm_year, poly3_GR_mm_year)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_GR_mm_year, gam_GR_mm_year, poly2_GR_mm_year, poly3_GR_mm_year, test = "F")

# Plot fitted values
plot(data_TOXEM$GR_mm_year, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$GR_mm_year, fitted(lm_GR_mm_year), col = "green", pch = 4)
points(data_TOXEM$GR_mm_year, fitted(gam_GR_mm_year), col = "blue", pch = 3)
points(data_TOXEM$GR_mm_year, fitted(poly2_GR_mm_year), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$GR_mm_year, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$GR_mm_year, predict(lm_GR_mm_year), col = "green")
points(data_TOXEM$GR_mm_year, predict(gam_GR_mm_year), col = "blue")
points(data_TOXEM$GR_mm_year, predict(poly2_GR_mm_year), col = "red")

# Check assumptions
plot(gam_GR_mm_year)
res <- resid(gam_GR_mm_year)
plot(fitted(gam_GR_mm_year), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
data_TOXEM <- na.omit(data_TOXEM[, c("GR_mm_year", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot(data_TOXEM)+ 
  aes(GR_mm_year, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_TOXEM$GR_mm_year, y=data_TOXEM$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_GR_mm_year <- lm(EROD_pmol_min_mg_prot ~ GR_mm_year + station * season * sex, data = data_TOXEM)
summary(lm_GR_mm_year)
# GR_mm_year significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_GR_mm_year <- gam(EROD_pmol_min_mg_prot ~ s(GR_mm_year) + station * season * sex, data = data_TOXEM)
summary(gam_GR_mm_year)
# GR_mm_year significantly affect the response in this model.

# Polynomial
poly2_GR_mm_year <- lm(EROD_pmol_min_mg_prot ~ poly(GR_mm_year, 2) + station * season * sex, data = data_TOXEM)
summary(poly2_GR_mm_year)
# GR_mm_year significantly affect the response in this model.

poly3_GR_mm_year <- lm(EROD_pmol_min_mg_prot ~ poly(GR_mm_year, 3) + station * season * sex, data = data_TOXEM)
summary(poly3_GR_mm_year)
# GR_mm_year does not significantly affect the response in this model.

### weight_g ***----------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = weight_g, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(weight_g ~ station + season + sex, data = data_TOXEM)
Anova(aov01, type = "III")

# Bar plot per station & season
ggplot(data_TOXEM) +
  aes(season, weight_g, fill = sex) +
  geom_boxplot() +
  facet_wrap(~station)

aov01 <- aov(weight_g ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in weight_g between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_TOXEM) +
  aes(weight_g, log_prop_hepatovac_area, col = season, group = season, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_TOXEM$weight_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_weight_g <- lm(log_prop_hepatovac_area ~ weight_g + station * season * sex, data = data_TOXEM)
summary(lm_weight_g)
# weight_g does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_weight_g <- gam(log_prop_hepatovac_area ~ s(weight_g) + station * season * sex, data = data_TOXEM)
summary(gam_weight_g)
plot(gam_weight_g, se = TRUE, scale = 0, pages = 1)
# weight_g does not significantly affect the response in this model.

# Polynomial
poly2_weight_g <- lm(log_prop_hepatovac_area ~ poly(weight_g, 2) + station * season * sex, data = data_TOXEM)
summary(poly2_weight_g)
# weight_g does not significantly affect the response in this model.

ggplot(data_TOXEM) +
  aes(weight_g, log_prop_hepatovac_area, col = season, group = season, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + 
  facet_wrap(~ station)

poly3_weight_g <- lm(log_prop_hepatovac_area ~ poly(weight_g, 3) + station * season * sex, data = data_TOXEM)
summary(poly3_weight_g)
# weight_g has a non-linear (cubic) significant effect on log_prop_hepatovac_area in this model.

# Model optimization
poly3_weight_g_b <- lm(log_prop_hepatovac_area ~ poly(weight_g, 3) + season, data = data_TOXEM)
summary(poly3_weight_g_b)
plot(poly3_weight_g_b)

residuals <- rstandard(poly3_weight_g_b)
outliers <- abs(residuals) > 3  # or 3 for stricter threshold
clean_data <- data_TOXEM[!outliers, ]

poly3_weight_g_c <- lm(log_prop_hepatovac_area ~ poly(weight_g, 3) + season, data = clean_data)
summary(poly3_weight_g_c)

ggplot(clean_data, aes(x = weight_g, y = log_prop_hepatovac_area, color = season, group = season)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, size = 1.5) +
  labs(title = "Polynomial Regression: log_prop_hepatovac_area ~ poly(weight_g, 3) + season", x = "Weight (g)", y = "log(prop_hepatovac_area)") +
  theme_minimal() + 
  facet_wrap(~ station)

# 2. Plot residuals vs fitted values for diagnostics
plot(poly3_weight_g_c, which = 1,
     main = "Residuals vs Fitted",
     sub.caption = "Checking for heteroscedasticity and outliers")

# 3. Plot residuals vs weight_g
ggplot(data.frame(weight_g = clean_data$weight_g,
                  residuals = residuals(poly3_weight_g_c),
                  season = clean_data$season),
       aes(x = weight_g, y = residuals, color = season)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Weight by Season",
       x = "Weight (g)", y = "Residuals") +
  theme_minimal()

# 4. Partial effect plot for poly(weight_g, 3)
library(effects)
plot(allEffects(poly3_weight_g_c), main = "Partial Effects Plot")

# Compare models
vif(lm_weight_g)
vif(gam_weight_g)
vif(poly2_weight_g)
vif(poly3_weight_g_b)
vif(poly3_weight_g_c)
# vif < 10 -> no collinearity between variables

summary(lm_weight_g)$r.squared
summary(gam_weight_g)$r.sq
summary(poly2_weight_g)$r.squared
summary(poly3_weight_g)$r.squared

AIC(lm_weight_g, gam_weight_g, poly2_weight_g, poly3_weight_g)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_weight_g, gam_weight_g, poly2_weight_g, poly3_weight_g, test = "F")

# Plot fitted values
plot(data_TOXEM$weight_g, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$weight_g, fitted(lm_weight_g), col = "green", pch = 4)
points(data_TOXEM$weight_g, fitted(gam_weight_g), col = "blue", pch = 3)
points(data_TOXEM$weight_g, fitted(poly2_weight_g), col = "red", pch = 2)
points(data_TOXEM$weight_g, fitted(poly3_weight_g), col = "orange", pch = 5)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$weight_g, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$weight_g, predict(lm_weight_g), col = "green")
points(data_TOXEM$weight_g, predict(gam_weight_g), col = "blue")
points(data_TOXEM$weight_g, predict(poly2_weight_g), col = "red")

# Check assumptions
plot(lm_weight_g)
res <- resid(lm_weight_g)
plot(fitted(lm_weight_g), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly3_weight_g model show a significant effect of weight_g on log_prop_hepatovac_area.
# poly3_weight_g

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
data_TOXEM <- na.omit(data_TOXEM[, c("weight_g", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot(data_TOXEM)+ 
  aes(weight_g, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_TOXEM$weight_g, y=data_TOXEM$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_weight_g <- lm(EROD_pmol_min_mg_prot ~ weight_g + station * season * sex, data = data_TOXEM)
summary(lm_weight_g)
# weight_g significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_weight_g <- gam(EROD_pmol_min_mg_prot ~ s(weight_g) + station * season * sex, data = data_TOXEM)
summary(gam_weight_g)
# weight_g significantly affect the response in this model.

# Polynomial
poly2_weight_g <- lm(EROD_pmol_min_mg_prot ~ poly(weight_g, 2) + station * season * sex, data = data_TOXEM)
summary(poly2_weight_g)
# weight_g does not significantly affect the response in this model.

poly3_weight_g <- lm(EROD_pmol_min_mg_prot ~ poly(weight_g, 3) + station * season * sex, data = data_TOXEM)
summary(poly3_weight_g)
# weight_g does not significantly affect the response in this model.

### carcasse_weight_g *** --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(carcasse_weight_g))
#data_clean <- na.omit(data_TOXEM[, c("carcasse_weight_g", "log_prop_hepatovac_area", "station", "season","sex")])

# Bar plot per station & season
ggplot(data_clean) +
  aes(season, carcasse_weight_g, fill = sex) +
  geom_col(stat = "identity", position = "dodge") +
  facet_wrap(~station)

aov01 <- aov(carcasse_weight_g ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in carcasse_weight_g between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_clean) +
  aes(carcasse_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_clean$carcasse_weight_g, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_carcasse_weight <- lm(log_prop_hepatovac_area ~ carcasse_weight_g + station * season * sex, data = data_clean)
summary(lm_carcasse_weight)
# carcasse_weight_g does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(carcasse_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_carcasse_weight <- gam(log_prop_hepatovac_area ~ s(carcasse_weight_g) + station * season * sex, data = data_clean)
summary(gam_carcasse_weight)
# carcasse_weight_g does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(carcasse_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_carcasse_weight <- lm(log_prop_hepatovac_area ~ poly(carcasse_weight_g, 2) + station * season * sex, data = data_clean)
summary(poly2_carcasse_weight)

# carcasse_weight_g does not significantly affect the response in this model.
ggplot(data_clean) +
  aes(carcasse_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_carcasse_weight <- lm(log_prop_hepatovac_area ~ poly(carcasse_weight_g, 3) + station * season * sex, data = data_clean)
summary(poly3_carcasse_weight)
# carcasse_weight_g does not significantly affect the response in this model.

# Compare models
vif(lm_carcasse_weight)
vif(gam_carcasse_weight)
vif(poly2_carcasse_weight)
vif(poly3_carcasse_weight)
# vif < 10 -> no collinearity between variables

summary(lm_carcasse_weight)$r.squared
summary(gam_carcasse_weight)$r.sq
summary(poly2_carcasse_weight)$r.squared
summary(poly3_carcasse_weight)$r.squared

AIC(lm_carcasse_weight, gam_carcasse_weight, poly2_carcasse_weight, poly3_carcasse_weight)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_carcasse_weight, gam_carcasse_weight, poly2_carcasse_weight, poly3_carcasse_weight, test = "F")
# Plot fitted values

plot(data_TOXEM$carcasse_weight_g, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$carcasse_weight_g, fitted(lm_carcasse_weight), col = "green", pch = 4)
points(data_clean$carcasse_weight_g, fitted(gam_carcasse_weight), col = "blue", pch = 3)
points(data_clean$carcasse_weight_g, fitted(poly2_carcasse_weight), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$carcasse_weight_g, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$carcasse_weight_g, predict(lm_carcasse_weight), col = "green")
points(data_TOXEM$carcasse_weight_g, predict(gam_carcasse_weight), col = "blue")
points(data_TOXEM$carcasse_weight_g, predict(poly2_carcasse_weight), col = "red")

# Check assumptions
plot(lm_carcasse_weight)
res <- resid(lm_carcasse_weight)
plot(fitted(lm_carcasse_weight), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### liver_weight_g ***--------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = liver_weight_g, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(liver_weight_g))
#data_clean <- na.omit(data_TOXEM[, c("liver_weight_g", "log_prop_hepatovac_area", "station", "season","sex")])

# Bar plot per station & season
ggplot(data_clean) +
  aes(season, liver_weight_g, fill = sex) +
  geom_col(stat = "identity", position = "dodge") +
  facet_wrap(~station)

aov01 <- aov(liver_weight_g ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in liver_weight_g between fishes from Seine and Canche and between females and males.

# Linear regression
ggplot(data_clean) +
  aes(liver_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_clean$liver_weight_g, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_liver_weight_g <- lm(log_prop_hepatovac_area ~ liver_weight_g + station * season * sex, data = data_clean)
summary(lm_liver_weight_g)
# liver_weight_g significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(liver_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_liver_weight_g <- gam(log_prop_hepatovac_area ~ s(liver_weight_g) + station * season * sex, data = data_clean)
summary(gam_liver_weight_g)
# liver_weight_g does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(liver_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_liver_weight_g <- lm(log_prop_hepatovac_area ~ poly(liver_weight_g, 2) + station * season * sex, data = data_clean)
summary(poly2_liver_weight_g)
# poly2_liver_weight_g model shows a significant affect of liver_weight_g on log_prop_hepatovac_area.

ggplot(data_clean) +
  aes(liver_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_liver_weight_g <- lm(log_prop_hepatovac_area ~ poly(liver_weight_g, 3) + station * season * sex, data = data_clean)
summary(poly3_liver_weight_g)
# poly3_liver_weight_g model shows a significant affect of liver_weight_g on log_prop_hepatovac_area.

# Compare models
vif(lm_liver_weight_g)
vif(gam_liver_weight_g)
vif(poly2_liver_weight_g)
vif(poly3_liver_weight_g)
# vif < 10 -> no collinearity between variables

summary(lm_liver_weight_g)$r.squared
summary(gam_liver_weight_g)$r.sq
summary(poly2_liver_weight_g)$r.squared
summary(poly3_liver_weight_g)$r.squared

AIC(lm_liver_weight_g, gam_liver_weight_g, poly2_liver_weight_g, poly3_liver_weight_g)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_liver_weight_g, gam_liver_weight_g, poly2_liver_weight_g, poly3_liver_weight_g, test = "F")

# Plot fitted values
plot(data_TOXEM$liver_weight_g, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$liver_weight_g, fitted(lm_liver_weight_g), col = "green", pch = 4)
points(data_clean$liver_weight_g, fitted(gam_liver_weight_g), col = "blue", pch = 3)
points(data_clean$liver_weight_g, fitted(poly3_liver_weight_g), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$liver_weight_g, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$liver_weight_g, predict(lm_liver_weight_g), col = "green")
points(data_TOXEM$liver_weight_g, predict(gam_liver_weight_g), col = "blue")
points(data_TOXEM$liver_weight_g, predict(poly2_liver_weight_g), col = "red")

# Check assumptions
plot(poly3_liver_weight_g)
res <- resid(poly3_liver_weight_g)
plot(fitted(poly3_liver_weight_g), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly3_liver_weight_g model show a significant effect of liver_weight_g on log_prop_hepatovac_area.
# poly3_liver_weight_g

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(liver_weight_g))
#data_clean <- na.omit(data_TOXEM[, c("liver_weight_g", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot(data_clean)+ 
  aes(liver_weight_g, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_clean$liver_weight_g, y=data_clean$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_liver_weight_g <- lm(EROD_pmol_min_mg_prot ~ liver_weight_g + station * season * sex, data = data_clean)
summary(lm_liver_weight_g)
# liver_weight_g significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_liver_weight_g <- gam(EROD_pmol_min_mg_prot ~ s(liver_weight_g) + station * season * sex, data = data_clean)
summary(gam_liver_weight_g)
# liver_weight_g significantly affect the response in this model.

# Polynomial
poly2_liver_weight_g <- lm(EROD_pmol_min_mg_prot ~ poly(liver_weight_g, 2) + station * season * sex, data = data_clean)
summary(poly2_liver_weight_g)
# liver_weight_g does not significantly affect the response in this model.

poly3_liver_weight_g <- lm(EROD_pmol_min_mg_prot ~ poly(liver_weight_g, 3) + station * season * sex, data = data_clean)
summary(poly3_liver_weight_g)
# liver_weight_g does not significantly affect the response in this model.

### HSI ***--------------------------------------------------------------
ggplot(data_TOXEM) +
  aes(x = season, y = HSI, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_TOXEM %>%
  group_by(station, season, sex) %>%
  summarise(mean_HSI = mean(HSI))

aov01 <- aov(HSI ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in HSI between fishes from Seine and Canche.

# Linear regression
(HSI <- ggplot(data_TOXEM) +
  aes(HSI, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station))

cor.test(x = data_TOXEM$HSI, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# statistically significant positive correlation

lm_HSI <- lm(log_prop_hepatovac_area ~ HSI, data = data_TOXEM)
summary(lm_HSI)

lm_HSI <- lm(log_prop_hepatovac_area ~ HSI + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_HSI)
# lm_HSI have a significant effect on log_prop_hepatovac_area in this model.

lm_HSI_simply <- lm(log_prop_hepatovac_area ~ HSI + season, data = data_TOXEM)
summary(lm_HSI_simply)

AIC(lm_HSI, lm_HSI_simply)

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(HSI, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_HSI <- gam(log_prop_hepatovac_area ~ s(HSI) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_HSI)
# gam_HSI shows a significant effect on log_prop_hepatovac_area in this model.

# Polynomial
ggplot(data_TOXEM) +
  aes(HSI, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_HSI <- lm(log_prop_hepatovac_area ~ poly(HSI, 2) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(poly2_HSI)
# poly2_HSI model shows a significant effect on log_prop_hepatovac_area in this model.

ggplot(data_TOXEM) +
  aes(HSI, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_HSI <- lm(log_prop_hepatovac_area ~ poly(HSI, 3) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(poly3_HSI)
# poly3_HSI model shows a significant effect on log_prop_hepatovac_area in this model.

AIC(lm_HSI, poly2_HSI, poly3_HSI)

# Compare models
vif(lm_HSI)
vif(gam_HSI)
vif(poly2_HSI)
vif(poly3_HSI)
# vif < 10 -> no collinearity between variables

summary(lm_HSI)$r.squared
summary(gam_HSI)$r.sq
summary(poly2_HSI)$r.squared
summary(poly3_HSI)$r.squared

AIC(lm_HSI, gam_HSI, poly2_HSI, poly3_HSI)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_HSI, gam_HSI, poly2_HSI, poly3_HSI, test = "F")

# Plot fitted values
plot(data_TOXEM$HSI, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$HSI, fitted(lm_HSI), col = "green", pch = 4)
points(data_TOXEM$HSI, fitted(gam_HSI), col = "blue", pch = 3)
points(data_TOXEM$HSI, fitted(poly2_HSI), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$HSI, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$HSI, predict(lm_HSI), col = "green")
points(data_TOXEM$HSI, predict(gam_HSI), col = "blue")
points(data_TOXEM$HSI, predict(poly2_HSI), col = "red")

# Check assumptions
plot(gam_HSI)
res <- resid(gam_HSI)
plot(fitted(gam_HSI), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# gam_HSI model show a significant effect of HSI on log_prop_hepatovac_area.
# gam_HSI

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
data_TOXEM <- na.omit(data_TOXEM[, c("HSI", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot(data_TOXEM)+ 
  aes(HSI, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_TOXEM$HSI, y=data_TOXEM$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_HSI <- lm(EROD_pmol_min_mg_prot ~ HSI + station * season * sex, data = data_TOXEM)
summary(lm_HSI)
# HSI significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_HSI <- gam(EROD_pmol_min_mg_prot ~ s(HSI) + station * season * sex, data = data_TOXEM)
summary(gam_HSI)
# HSI significantly affect the response in this model.

# Polynomial
poly2_HSI <- lm(EROD_pmol_min_mg_prot ~ poly(HSI, 2) + station * season * sex, data = data_TOXEM)
summary(poly2_HSI)
# HSI significantly affect the response in this model.

poly3_HSI <- lm(EROD_pmol_min_mg_prot ~ poly(HSI, 3) + station * season * sex, data = data_TOXEM)
summary(poly3_HSI)
# HSI does not significantly affect the response in this model.

# Including HSI as a covariate accounts for the significant difference in mean HSI between the two fish populations.
# HSI will be used as a covariate in every model 

### TG_FS_ratio --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(TG_FS_ratio))

ggplot(data_TOXEM) +
  aes(x = season, y = TG_FS_ratio, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 20)) +
  #geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_TG_FS_ratio = mean(TG_FS_ratio),
            sd_TG_FS_ratio = sd(TG_FS_ratio))

aov01 <- aov(TG_FS_ratio ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in TG_FS_ratio between fishes from Seine and Canche, and between females and males.


# Linear regression
ggplot(data_TOXEM) +
  aes(TG_FS_ratio, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$TG_FS_ratio, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_TG_FS <- lm(log_prop_hepatovac_area ~ TG_FS_ratio + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_TG_FS)
# TG_FS_ratio does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(TG_FS_ratio, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_TG_FS <- gam(log_prop_hepatovac_area ~ s(TG_FS_ratio) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_TG_FS)
# TG_FS_ratio does not significantly affect the response in this model.

# Polynomial
data_clean <- subset(data_TOXEM, !is.na(TG_FS_ratio))
#data_clean <- na.omit(data_TOXEM[, c("TG_FS_ratio", "log_prop_hepatovac_area", "station", "season", "sex", "age_month")]) 

ggplot(data_TOXEM) +
  aes(TG_FS_ratio, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_lipids_TG_FS <- lm(log_prop_hepatovac_area ~ poly(TG_FS_ratio, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_TG_FS)
# TG_FS_ratio does not significantly affect the response in this model.

ggplot(data_TOXEM) +
  aes(TG_FS_ratio, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_lipids_TG_FS <- lm(log_prop_hepatovac_area ~ poly(TG_FS_ratio, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_TG_FS)
# TG_FS_ratio does not significantly affect the response in this model.

# Compare models
vif(lm_lipids_TG_FS)
vif(gam_lipids_TG_FS)
vif(poly2_lipids_TG_FS)
vif(poly3_lipids_TG_FS)
# vif < 10 -> no collinearity between variables

summary(lm_lipids_TG_FS)$r.squared
summary(gam_lipids_TG_FS)$r.sq
summary(poly2_lipids_TG_FS)$r.squared
summary(poly3_lipids_TG_FS)$r.squared

AIC(lm_lipids_TG_FS, gam_lipids_TG_FS, poly2_lipids_TG_FS, poly3_lipids_TG_FS)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_lipids_TG_FS, gam_lipids_TG_FS, poly2_lipids_TG_FS, poly3_lipids_TG_FS, test = "F")
# Plot fitted values
plot(data_TOXEM$TG_FS_ratio, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$TG_FS_ratio, fitted(lm_lipids_TG_FS), col = "green", pch = 4)
points(data_TOXEM$TG_FS_ratio, fitted(gam_lipids_TG_FS), col = "blue", pch = 3)
points(data_TOXEM$TG_FS_ratio, fitted(poly2_lipids_TG_FS), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))
plot(data_TOXEM$TG_FS_ratio, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$TG_FS_ratio, predict(lm_lipids_TG_FS), col = "green")
points(data_TOXEM$TG_FS_ratio, predict(gam_lipids_TG_FS), col = "blue")
points(data_TOXEM$TG_FS_ratio, predict(poly2_lipids_TG_FS), col = "red")
# Check assumptions
plot(lm_lipids_TG_FS)
res <- resid(lm_lipids_TG_FS)
plot(fitted(lm_lipids_TG_FS), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# gam_lipids_TG_FS


## plot_fitness_data_TOXEM -------------------------------------------------
ggarrange(K_index, HSI, age_month, ncol = 2, nrow = 2, labels = c("A","B", "C"))

## Reproductive_marker -----------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

data_female <- data_TOXEM %>% 
  filter(sex == 'F')

data_male <- data_TOXEM %>% 
  filter(sex == 'M')

### GSI_females --------------------------------------------------------------
ggplot(data_female) +
  aes(x = weight_g, y = gonad_weight_g, col = season) +
  geom_point() +
  facet_wrap(~ station)

ggplot(data_female) +
  aes(x = season, y = GSI) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_female %>% 
  group_by(season) %>% 
  summarise(mean_GSI = mean(GSI), 
            sd_GSI = sd(GSI))

aov01 <- aov(GSI ~ station * season, data = data_female)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in GSI between female fish from Seine and Canche and between female fish sampled in Winter vs Summer.

# Linear regression
ggplot(data_female) +
  aes(GSI, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_female$GSI, y = data_female$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no significant correlation

lm_GSI <- lm(log_prop_hepatovac_area ~ GSI + station * season, data = data_female)
summary(lm_GSI)
# GSI does not significantly affect the response in this model.

# Non linear regression
# GAM
gam_GSI <- gam(log_prop_hepatovac_area ~ s(GSI) + station * season, data = data_female)
summary(gam_GSI)
# GSI does not significantly affect the response in this model.

poly2_GSI <- lm(log_prop_hepatovac_area ~ poly(GSI, 2) + station * season, data = data_female)
summary(poly2_GSI)
# GSI does not significantly affect the response in this model.

poly3_GSI <- lm(log_prop_hepatovac_area ~ poly(GSI, 3) + station * season, data = data_female)
summary(poly3_GSI)
# GSI does not significantly affect the response in this model.

### gonad_weight_g_females *** --------------------------------------------------------------
ggplot(data_female) +
  aes(x = season, y = gonad_weight_g) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

# Remove rows with missing values in relevant columns
data_clean <- subset(data_female, !is.na(gonad_weight_g))
#data_clean <- na.omit(data_female[, c("gonad_weight_g", "log_prop_hepatovac_area", "station", "season","sex")])

# Bar plot per station & season
ggplot(data_clean) +
  aes(season, gonad_weight_g, fill = sex) +
  geom_col(stat = "identity", position = "dodge") +
  facet_wrap(~station)

aov01 <- aov(gonad_weight_g ~ station * season * sex, data = data_female)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in gonad_weight_g between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_clean$gonad_weight_g, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no significant correlation

lm_gonad_weight_g <- lm(log_prop_hepatovac_area ~ gonad_weight_g + station * season, data = data_clean)
summary(lm_gonad_weight_g)
# gonad_weight_g significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_gonad_weight_g <- gam(log_prop_hepatovac_area ~ s(gonad_weight_g) + station * season * sex, data = data_clean)
summary(gam_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_gonad_weight_g <- lm(log_prop_hepatovac_area ~ poly(gonad_weight_g, 2) + station * season, data = data_clean)
summary(poly2_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_gonad_weight_g <- lm(log_prop_hepatovac_area ~ poly(gonad_weight_g, 3) + station * season, data = data_clean)
summary(poly3_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

# Compare models
vif(lm_gonad_weight_g)
vif(gam_gonad_weight_g)
vif(poly2_gonad_weight_g)
vif(poly3_gonad_weight_g)
# vif < 10 -> no collinearity between variables

summary(lm_gonad_weight_g)$r.squared
summary(gam_gonad_weight_g)$r.sq
summary(poly2_gonad_weight_g)$r.squared
summary(poly3_gonad_weight_g)$r.squared

AIC(lm_gonad_weight_g, gam_gonad_weight_g, poly2_gonad_weight_g, poly3_gonad_weight_g)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_gonad_weight_g, gam_gonad_weight_g, poly2_gonad_weight_g, poly3_gonad_weight_g, test = "F")

# Plot fitted values
plot(data_female$gonad_weight_g, data_female$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$gonad_weight_g, fitted(lm_gonad_weight_g), col = "green", pch = 4)
points(data_clean$gonad_weight_g, fitted(gam_gonad_weight_g), col = "blue", pch = 3)
points(data_clean$gonad_weight_g, fitted(poly2_gonad_weight_g), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_female$gonad_weight_g, data_female$log_prop_hepatovac_area)
points(data_female$gonad_weight_g, predict(lm_gonad_weight_g), col = "green")
points(data_female$gonad_weight_g, predict(gam_gonad_weight_g), col = "blue")
points(data_female$gonad_weight_g, predict(poly2_gonad_weight_g), col = "red")

# Check assumptions
plot(lm_gonad_weight_g)
res <- resid(lm_gonad_weight_g)
plot(fitted(lm_gonad_weight_g), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### sexual_maturation_total_female --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_female, !is.na(sexual_maturation_total))

ggplot(data_female) +
  aes(x = season, y = sexual_maturation_total) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

#kruskal.test(sexual_maturation_total ~ station, data = data_female)

aov01 <- aov(sexual_maturation_total ~ station * season, data = data_female)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in sexual_maturation_total between fishes from the different station and season.

# Linear regression
ggplot(data_female) +
  aes(sexual_maturation_total, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_female$sexual_maturation_total, y = data_female$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# not significant 

lm_sexual_maturation_total <- lm(log_prop_hepatovac_area ~ sexual_maturation_total + station * season + poly(age_month, 2), data = data_female)
summary(lm_sexual_maturation_total)
# not significant

# Non linear regression
# GAM
ggplot(data_female) +
  aes(sexual_maturation_total, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_sexual_maturation_total <- gam(log_prop_hepatovac_area ~ s(sexual_maturation_total) + station * season + poly(age_month, 2), data = data_female)
summary(gam_sexual_maturation_total)
# not significant

# Polynomial
ggplot(data_female) +
  aes(sexual_maturation_total, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_sexual_maturation_total <- lm(log_prop_hepatovac_area ~ poly(sexual_maturation_total, 2) + season, data = data_female)
summary(poly2_sexual_maturation_total)
# significant

ggplot(data_female) +
  aes(sexual_maturation_total, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_sexual_maturation_total <- lm(log_prop_hepatovac_area ~ poly(sexual_maturation_total, 3) + poly(age_month, 2) + season, data = data_female)
summary(poly3_sexual_maturation_total)
# significant

vif(poly3_sexual_maturation_total)

# Compare models
vif(lm_sexual_maturation_total)
vif(gam_sexual_maturation_total)
vif(poly2_sexual_maturation_total)
vif(poly3_sexual_maturation_total)
# vif < 10 -> no collinearity between variables

summary(lm_sexual_maturation_total)$r.squared
summary(gam_sexual_maturation_total)$r.sq
summary(poly2_sexual_maturation_total)$r.squared
summary(poly3_sexual_maturation_total)$r.squared

AIC(lm_sexual_maturation_total, gam_sexual_maturation_total, poly2_sexual_maturation_total, poly3_sexual_maturation_total)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_sexual_maturation_total, gam_sexual_maturation_total, poly2_sexual_maturation_total, poly3_sexual_maturation_total, test = "F")

# Plot fitted values
plot(data_female$sexual_maturation_total, data_female$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_female$sexual_maturation_total, fitted(lm_sexual_maturation_total), col = "green", pch = 4)
points(data_female$sexual_maturation_total, fitted(gam_sexual_maturation_total), col = "blue", pch = 3)
points(data_female$sexual_maturation_total, fitted(poly2_sexual_maturation_total), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_female$sexual_maturation_total, data_female$log_prop_hepatovac_area)
points(data_female$sexual_maturation_total, predict(lm_sexual_maturation_total), col = "green")
points(data_female$sexual_maturation_total, predict(gam_sexual_maturation_total), col = "blue")
points(data_female$sexual_maturation_total, predict(poly2_sexual_maturation_total), col = "red")

# Check assumptions
plot(gam_sexual_maturation_total)
res <- resid(gam_sexual_maturation_total)
plot(fitted(gam_sexual_maturation_total), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# gam_sexual_maturation_total model show a significant effect of sexual_maturation_total on log_prop_hepatovac_area.
# gam_sexual_maturation_total

# Comparison with EROD_pmol_min_mg_prot
# Remove rows with missing values in relevant columns
data_female <- na.omit(data_female[, c("sexual_maturation_total", "EROD_pmol_min_mg_prot", "station", "season","sex")])

# Linear regression
ggplot(data_female)+ 
  aes(sexual_maturation_total, EROD_pmol_min_mg_prot, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x=data_female$sexual_maturation_total, y=data_female$EROD_pmol_min_mg_prot, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_sexual_maturation_total <- lm(EROD_pmol_min_mg_prot ~ sexual_maturation_total + station * season * sex, data = data_female)
summary(lm_sexual_maturation_total)
# sexual_maturation_total significantly affect the response but the model is not significant as a whole

# Non linear regression
# GAM
gam_sexual_maturation_total <- gam(EROD_pmol_min_mg_prot ~ s(sexual_maturation_total) + station * season * sex, data = data_female)
summary(gam_sexual_maturation_total)
# sexual_maturation_total significantly affect the response in this model.

# Polynomial
poly2_sexual_maturation_total <- lm(EROD_pmol_min_mg_prot ~ poly(sexual_maturation_total, 2) + station * season * sex, data = data_female)
summary(poly2_sexual_maturation_total)
# sexual_maturation_total significantly affect the response in this model.

poly3_sexual_maturation_total <- lm(EROD_pmol_min_mg_prot ~ poly(sexual_maturation_total, 3) + station * season * sex, data = data_female)
summary(poly3_sexual_maturation_total)
# sexual_maturation_total does not significantly affect the response in this model.

### E2 (oestradiol)_female *** --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_female, !is.na(E2))
# Remove rows with E2 values for sex == M -> measured only in females
#data_clean <- data_clean[!(data_clean$sex == "M" & !is.na(data_clean$E2)), ]

ggplot(data_female) +
  aes(x = season, y = E2) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>% 
  group_by(season) %>% 
  summarise(mean_E2 = mean(E2), 
            sd_E2 = sd(E2))

aov01 <- aov(E2 ~ station * season, data = data_female)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in E2 level between fishes from winter and summer. No significant difference between stations.

# Linear regression
(E2 <- ggplot(data_female) +
  aes(E2, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station))

cor.test(x = data_female$E2, y = data_female$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# negative statistically significant correlation

lm_E2 <- lm(log_prop_hepatovac_area ~ E2, data = data_female)
summary(lm_E2)
# E2 does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_female) +
  aes(E2, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_E2 <- gam(log_prop_hepatovac_area ~ s(E2) + season + poly(age_month, 2), data = data_female)
summary(gam_E2)
# E2 does not significantly affect the response in this model.

# Polynomial
ggplot(data_female) +
  aes(E2, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_E2 <- lm(log_prop_hepatovac_area ~ poly(E2, 2) + season + poly(age_month, 2), data = data_clean)
summary(poly2_E2)
# E2 does not significantly affect the response in this model.

ggplot(data_female) +
  aes(E2, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_E2 <- lm(log_prop_hepatovac_area ~ poly(E2, 3) + season + poly(age_month, 2), data = data_clean)
summary(poly3_E2)
# E2 does not significantly affect the response in this model.

# Compare models
vif(lm_E2)
vif(gam_E2)
vif(poly2_E2)
vif(poly3_E2)
# vif < 10 -> no collinearity between variables

summary(lm_E2)$r.squared
summary(gam_E2)$r.sq
summary(poly2_E2)$r.squared
summary(poly3_E2)$r.squared

AIC(lm_E2, gam_E2, poly2_E2, poly3_E2)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_E2, gam_E2, poly2_E2, poly3_E2, test = "F")

# Plot fitted values
plot(data_TOXEM$E2, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_female$E2, fitted(lm_E2), col = "green", pch = 4)
points(data_female$E2, fitted(gam_E2), col = "blue", pch = 3)
points(data_female$E2, fitted(poly2_E2), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$E2, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$E2, predict(lm_E2), col = "green")
points(data_TOXEM$E2, predict(gam_E2), col = "blue")
points(data_TOXEM$E2, predict(poly2_E2), col = "red")

# Check assumptions
plot(lm_E2)
res <- resid(lm_E2)
plot(fitted(lm_E2), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### GSI_male *** --------------------------------------------------------------
ggplot(data_male) +
  aes(x = weight_g, y = gonad_weight_g, col = season) +
  geom_point() +
  facet_wrap(~ station)

ggplot(data_male) +
  aes(x = season, y = GSI) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_male %>%
  group_by(station) %>%
  summarise(mean_GSI = mean(GSI),
            sd_GSI = sd(GSI))

aov01 <- aov(GSI ~ station * season, data = data_male)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in GSI between female fish from Seine and Canche and between female fish sampled in Winter vs Summer.

# Linear regression
(GSI_male <- ggplot(data_male) +
  aes(GSI, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station))

cor.test(x = data_male$GSI, y = data_male$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# significant negative correlation

lm_GSI <- lm(log_prop_hepatovac_area ~ GSI, data = data_male)
summary(lm_GSI)
# GSI does not significantly affect the response in this model.

# Non linear regression
# GAM
gam_GSI <- gam(log_prop_hepatovac_area ~ s(GSI) + station * season, data = data_male)
summary(gam_GSI)
# GSI does not significantly affect the response in this model.

poly2_GSI <- lm(log_prop_hepatovac_area ~ poly(GSI, 2) + station * season, data = data_male)
summary(poly2_GSI)
# GSI does not significantly affect the response in this model.

poly3_GSI <- lm(log_prop_hepatovac_area ~ poly(GSI, 3) + station * season, data = data_male)
summary(poly3_GSI)
# GSI does not significantly affect the response in this model.

### gonad_weight_g_males *** --------------------------------------------------------------
ggplot(data_male) +
  aes(x = season, y = gonad_weight_g) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(gonad_weight_g))
#data_clean <- na.omit(data_male[, c("gonad_weight_g", "log_prop_hepatovac_area", "station", "season","sex")])

aov01 <- aov(gonad_weight_g ~ station * season * sex, data = data_male)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in gonad_weight_g between fishes from Seine and Canche, between fishes sampled in Winter vs Summer and between females and males.

# Linear regression
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_clean$gonad_weight_g, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no significant correlation

lm_gonad_weight_g <- lm(log_prop_hepatovac_area ~ gonad_weight_g + station * season, data = data_clean)
summary(lm_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_gonad_weight_g <- gam(log_prop_hepatovac_area ~ s(gonad_weight_g) + station * season, data = data_clean)
summary(gam_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_gonad_weight_g <- lm(log_prop_hepatovac_area ~ poly(gonad_weight_g, 2) + station * season, data = data_clean)
summary(poly2_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(gonad_weight_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_gonad_weight_g <- lm(log_prop_hepatovac_area ~ poly(gonad_weight_g, 3) + station * season, data = data_clean)
summary(poly3_gonad_weight_g)
# gonad_weight_g does not significantly affect the response in this model.

# Compare models
vif(lm_gonad_weight_g)
vif(gam_gonad_weight_g)
vif(poly2_gonad_weight_g)
vif(poly3_gonad_weight_g)
# vif < 10 -> no collinearity between variables

summary(lm_gonad_weight_g)$r.squared
summary(gam_gonad_weight_g)$r.sq
summary(poly2_gonad_weight_g)$r.squared
summary(poly3_gonad_weight_g)$r.squared

AIC(lm_gonad_weight_g, gam_gonad_weight_g, poly2_gonad_weight_g, poly3_gonad_weight_g)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_gonad_weight_g, gam_gonad_weight_g, poly2_gonad_weight_g, poly3_gonad_weight_g, test = "F")

# Plot fitted values
plot(data_male$gonad_weight_g, data_male$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$gonad_weight_g, fitted(lm_gonad_weight_g), col = "green", pch = 4)
points(data_clean$gonad_weight_g, fitted(gam_gonad_weight_g), col = "blue", pch = 3)
points(data_clean$gonad_weight_g, fitted(poly2_gonad_weight_g), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_male$gonad_weight_g, data_male$log_prop_hepatovac_area)
points(data_male$gonad_weight_g, predict(lm_gonad_weight_g), col = "green")
points(data_male$gonad_weight_g, predict(gam_gonad_weight_g), col = "blue")
points(data_male$gonad_weight_g, predict(poly2_gonad_weight_g), col = "red")

# Check assumptions
plot(lm_gonad_weight_g)
res <- resid(lm_gonad_weight_g)
plot(fitted(lm_gonad_weight_g), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### COMET_OTM_%_tail_intensity_male ---------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_male, !is.na(COMET_OTM))

ggplot(data_male) +
  aes(x = season, y = COMET_OTM) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

# Linear regression
ggplot(data_male) +
  aes(COMET_OTM, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = F)

cor.test(x = data_male$COMET_OTM, y = data_male$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_COMET_OTM <- lm(log_prop_hepatovac_area ~ COMET_OTM + poly(age_month, 2), data = data_male)
summary(lm_COMET_OTM)
# COMET_OTM does not significantly affect the response in this model.

# Non linear regression
gam_COMET_OTM <- gam(log_prop_hepatovac_area ~ s(COMET_OTM) + poly(age_month, 2), data = data_male)
summary(gam_COMET_OTM)
# COMET_OTM does not significantly affect the response in this model.

# Polynomial
poly2_COMET_OTM <- lm(log_prop_hepatovac_area ~ poly(COMET_OTM, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_COMET_OTM)
# COMET_OTM does not significantly affect the response in this model.

poly3_COMET_OTM <- lm(log_prop_hepatovac_area ~ poly(COMET_OTM, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_COMET_OTM)
# COMET_OTM does not significantly affect the response in this model.

### eleven_KT_male --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_male, !is.na(eleven_KT))

ggplot(data_male) +
  aes(x = season, y = eleven_KT) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_eleven_KT = mean(eleven_KT),
            sd_eleven_KT = sd(eleven_KT))

aov01 <- aov(eleven_KT ~ station * season, data = data_male)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in eleven_KT level between fishes from winter and summer. No significant difference between stations.

# Linear regression
ggplot(data_clean) +
  aes(eleven_KT, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_clean$eleven_KT, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_eleven_KT <- lm(log_prop_hepatovac_area ~ eleven_KT + station * season + poly(age_month, 2), data = data_clean)
summary(lm_eleven_KT)
# eleven_KT does not significantly affect the response in this model.

# Non linear regression
# GAM
gam_eleven_KT <- gam(log_prop_hepatovac_area ~ s(eleven_KT) + station * season + poly(age_month, 2), data = data_clean)
summary(gam_eleven_KT)
# eleven_KT does not significantly affect the response in this model.

# Polynomial
poly2_eleven_KT <- lm(log_prop_hepatovac_area ~ poly(eleven_KT, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_eleven_KT)
# eleven_KT does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(eleven_KT, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_eleven_KT <- lm(log_prop_hepatovac_area ~ poly(eleven_KT, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_eleven_KT)
# eleven_KT does not significantly affect the response in this model.

# Compare models
vif(lm_eleven_KT)
vif(gam_eleven_KT)
vif(poly2_eleven_KT)
vif(poly3_eleven_KT)
# vif < 10 -> no collinearity between variables

summary(lm_eleven_KT)$r.squared
summary(gam_eleven_KT)$r.sq
summary(poly2_eleven_KT)$r.squared
summary(poly3_eleven_KT)$r.squared

AIC(lm_eleven_KT, gam_eleven_KT, poly2_eleven_KT, poly3_eleven_KT)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_eleven_KT, gam_eleven_KT, poly2_eleven_KT, poly3_eleven_KT, test = "F")

# Plot fitted values
plot(data_TOXEM$eleven_KT, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$eleven_KT, fitted(lm_eleven_KT), col = "green", pch = 4)
points(data_clean$eleven_KT, fitted(gam_eleven_KT), col = "blue", pch = 3)
points(data_clean$eleven_KT, fitted(poly2_eleven_KT), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$eleven_KT, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$eleven_KT, predict(lm_eleven_KT), col = "green")
points(data_TOXEM$eleven_KT, predict(gam_eleven_KT), col = "blue")
points(data_TOXEM$eleven_KT, predict(poly2_eleven_KT), col = "red")

# Check assumptions
plot(lm_eleven_KT)
res <- resid(lm_eleven_KT)
plot(fitted(lm_eleven_KT), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### VTG_ng_ml_male --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_male, !is.na(VTG_ng_ml))

ggplot(data_male) +
  aes(x = season, y = VTG_ng_ml) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(VTG_ng_ml ~ station * season, data = data_male)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# no significant difference in VTG_ng_ml level between fishes from winter and summer or between stations.

# Linear regression
ggplot(data_male) +
  aes(VTG_ng_ml, log_prop_hepatovac_area, col = season, group = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_male$VTG_ng_ml, y = data_male$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_VTG <- lm(log_prop_hepatovac_area ~ VTG_ng_ml + station * season + poly(age_month, 2), data = data_clean)
summary(lm_VTG)
# VTG_ng_ml does not significantly affect the response in this model.

# Non linear regression
# GAM
gam_VTG <- gam(log_prop_hepatovac_area ~ s(VTG_ng_ml) + station * season + poly(age_month, 2), data = data_clean)
summary(gam_VTG)
# VTG_ng_ml does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(VTG_ng_ml, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_VTG <- lm(log_prop_hepatovac_area ~ poly(VTG_ng_ml, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_VTG)
# VTG_ng_ml does not significantly affect the response in this model.

poly3_VTG <- lm(log_prop_hepatovac_area ~ poly(VTG_ng_ml, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_VTG)
# VTG_ng_ml does not significantly affect the response in this model.

## plot_reproductive_marker -------------------------------------------------
ggarrange(E2, GSI_male, ncol = 2, nrow = 1, labels = c("A","B"))

## liver_histopathology ----------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

#data_TOXEM_histo_bis <- data_TOXEM %>%
#   select(station, season, sex, age_month, log_prop_hepatovac_area, non_specific_lesions_1, non_specific_lesions_2, non_specific_lesions_3, early_toxicopathic_non_neoplastic_lesions, FCA, benign_neoplasms, malignant_neoplasms)
 
#write_csv(data_TOXEM_histo_bis, "/home/valentin/Desktop/hepatovac/results/data_TOXEM_histo_bis.csv")

# Modify data_TOXEM_histo_bis externaly to have only one column for non_specific_lesions

data_histo <- read_csv("/home/valentin/Desktop/hepatovac/results/data_TOXEM_histo_bis.csv")

data_histo_longer <- data_histo %>% pivot_longer(cols=c("non_specific_lesions", "early_toxicopathic_non_neoplastic_lesions", "FCA", "benign_neoplasms", "malignant_neoplasms"),
                                        names_to='categories',
                                        values_to='lesions')
# Recode row with NAs 
data_histo_longer$categories[is.na(data_histo_longer$lesions)] <- NA
data_histo_longer <- data_histo_longer %>%
  group_by(log_prop_hepatovac_area) %>%
  filter(if (any(!is.na(lesions))) {!is.na(lesions)} else {row_number() == 1}) %>%
  ungroup()

counts <- data_histo_longer %>%
  mutate(categories = fct_explicit_na(categories, na_level = "NA")) %>%
  count(station, season, sex, categories)

(plot_categories <- ggplot(counts, aes(x = season, y = n, fill = categories)) +
  geom_col() +
  labs(x = "", y = "Count", fill = "Histological categories") +
  theme(text=element_text(size=20)) +
  facet_wrap(station~sex))

counts <- data_histo_longer %>%
  mutate(lesions = fct_explicit_na(lesions, na_level = "NA")) %>%
  count(station, season, sex, lesions)

(plot_lesions <- ggplot(counts, aes(x = season, y = n, fill = lesions)) +
  geom_col() +
  labs(x = "", y = "Count", fill = "Histological lesions") +
  theme(text=element_text(size=20)) +
  facet_wrap(station~sex))

ggarrange(plot_lesions, plot_categories, ncol = 2, nrow = 1, labels = c("A","B"))

### categories --------------------------------------------------------------
# Test if the prevalence of each categories is different between stations, seasons and sex
data_histo_categories <- data_histo %>%
  group_by(station, season, sex) %>%
  summarise(non_specific_lesions = sum(!is.na(non_specific_lesions)),
            early_toxicopathic_non_neoplastic_lesions = sum(!is.na(early_toxicopathic_non_neoplastic_lesions)),
            FCA = sum(!is.na(FCA)),
            benign_neoplasms = sum(!is.na(benign_neoplasms)),
            malignant_neoplasms = sum(!is.na(malignant_neoplasms)))

data_histo_categories <- merge(data_histo_categories, table_station_season_sex, by.x=c('station', 'season', 'sex'), by.y=c('station', 'season', 'sex'))
data_histo_categories <- data_histo_categories %>% relocate(count_samples, .after = sex)

data_histo_categories_prop <- data_histo_categories %>%
  mutate(prop_non_specific_lesions = non_specific_lesions*100/count_samples,
         prop_early_toxicopathic_non_neoplastic_lesions = early_toxicopathic_non_neoplastic_lesions*100/count_samples,
         prop_FCA = FCA*100/count_samples,
         prop_benign_neoplasms = benign_neoplasms*100/count_samples,
         prop_malignant_neoplasms = malignant_neoplasms*100/count_samples) %>%
  select(-c(count_samples, non_specific_lesions, early_toxicopathic_non_neoplastic_lesions, FCA, benign_neoplasms, malignant_neoplasms))

write_csv(data_histo_categories_prop, "/home/valentin/Desktop/hepatovac/results/data_histo_categories_prop.csv")

ggplot(data_histo_categories_prop) +
  aes(x = season, y = prop_non_specific_lesions, color = sex) +
  geom_point() +
  facet_wrap(~station)

ggplot(data_histo_categories_prop) +
  aes(x = season, y = prop_early_toxicopathic_non_neoplastic_lesions, color = sex) +
  geom_point() +
  facet_wrap(~station)

ggplot(data_histo_categories_prop) +
  aes(x = season, y = prop_FCA, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_non_specific_lesions ~ station, data = data_histo_categories_prop) # Significant difference detected in prevalence of non-specific lesions between stations, with higher prevalence in individuals from Seine compared to Canche.
kruskal.test(prop_non_specific_lesions ~ season, data = data_histo_categories_prop)
kruskal.test(prop_non_specific_lesions ~ sex, data = data_histo_categories_prop) 
dunnTest(prop_non_specific_lesions ~ station, data = data_histo_categories_prop, method="bonferroni")
ggplot(data_histo_categories_prop) +
  aes(x = season, y = prop_non_specific_lesions, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_early_toxicopathic_non_neoplastic_lesions ~ station, data = data_histo_categories_prop) 
kruskal.test(prop_early_toxicopathic_non_neoplastic_lesions ~ season, data = data_histo_categories_prop)
kruskal.test(prop_early_toxicopathic_non_neoplastic_lesions ~ sex, data = data_histo_categories_prop) 

kruskal.test(prop_FCA ~ station, data = data_histo_categories_prop) 
kruskal.test(prop_FCA ~ season, data = data_histo_categories_prop)
kruskal.test(prop_FCA ~ sex, data = data_histo_categories_prop) 

# Test if there is a relationship between the prevalence of each categories and log_prop_hepatovac_area 
data_histo_categories_prevalence <- data_histo
data_histo_categories_prevalence$prevalence_non_specific_lesions <- ifelse(!is.na(data_histo_categories_prevalence$non_specific_lesions), 1, 0)
data_histo_categories_prevalence$prevalence_early_toxicopathic_non_neoplastic_lesions <- ifelse(!is.na(data_histo_categories_prevalence$early_toxicopathic_non_neoplastic_lesions), 1, 0)
data_histo_categories_prevalence$prevalence_FCA <- ifelse(!is.na(data_histo_categories_prevalence$FCA), 1, 0)
data_histo_categories_prevalence$prevalence_benign_neoplasms <- ifelse(!is.na(data_histo_categories_prevalence$benign_neoplasms), 1, 0)
data_histo_categories_prevalence$prevalence_malignant_neoplasms <- ifelse(!is.na(data_histo_categories_prevalence$malignant_neoplasms), 1, 0)

priors <- get_prior(log_prop_hepatovac_area ~ prevalence_non_specific_lesions + 
                      prevalence_early_toxicopathic_non_neoplastic_lesions +
                      prevalence_FCA + 
                      #age_month + 
                      sex +
                      season + 
                      station,
                    data = data_histo_categories_prevalence, family = gaussian())

brm_histo_categories_model <- brms::brm(log_prop_hepatovac_area ~ prevalence_non_specific_lesions + 
                               prevalence_early_toxicopathic_non_neoplastic_lesions +
                               prevalence_FCA + 
                              #age_month + 
                               sex +
                               season + 
                               station, 
                             data=data_histo_categories_prevalence, 
                             prior = priors)

summary(brm_histo_categories_model)

### lesions -----------------------------------------------------------------
# Test if the prevalence of each lesions is different between stations, seasons and sex
data_histo_longer <- data_histo_longer %>% mutate(row_id = row_number())
data_histo_wider <- data_histo_longer %>% select(-categories)
data_histo_wider <- data_histo_wider %>% mutate(row_id = row_number())
data_histo_wider <- data_histo_wider %>%
  pivot_wider(names_from = lesions, values_from = lesions, id_cols = row_id)
data_histo_wider <- merge(data_histo_longer, data_histo_wider, by.x=c('row_id'), by.y=c('row_id'))
data_histo_wider <- data_histo_wider %>% select(-c(categories, lesions, row_id))
data_histo_wider <- data_histo_wider %>% select(-7)

# merge rows with duplicated log_prop_hepatovac_area values
data_histo_wider <- data_histo_wider %>%
  group_by(log_prop_hepatovac_area) %>%
  summarise(
    station = paste(unique(na.omit(station)), collapse = ";"),
    season = paste(unique(na.omit(season)), collapse = ";"),
    sex = paste(unique(na.omit(sex)), collapse = ";"),
    age_month = paste(unique(na.omit(age_month)), collapse = ";"),
    MMA = paste(unique(na.omit(MMA)), collapse = ";"),
    `Hydropic vacuolation` = paste(unique(na.omit(`Hydropic vacuolation`)), collapse = ";"),
    Lipidosis = paste(unique(na.omit(Lipidosis)), collapse = ";"),
    `Basophilic foci` = paste(unique(na.omit(`Basophilic foci`)), collapse = ";"),
    `Lymphocytic infiltration` = paste(unique(na.omit(`Lymphocytic infiltration`)), collapse = ";"),
    `Granuloma & Abscess` = paste(unique(na.omit(`Granuloma & Abscess`)), collapse = ";"),
    `Clear cell foci` = paste(unique(na.omit(`Clear cell foci`)), collapse = ";"),
    Necrosis = paste(unique(na.omit(Necrosis)), collapse = ";"),
    `Eosinophilic foci` = paste(unique(na.omit(`Eosinophilic foci`)), collapse = ";"),
    Phospholipidosis = paste(unique(na.omit(Phospholipidosis)), collapse = ";"),
    .groups = "drop"
  )

colnames(data_histo_wider)[colnames(data_histo_wider) == 'MMA'] <- 'MMA'
data_histo_wider$MMA[data_histo_wider$MMA == "MMA"] <- 1
data_histo_wider$MMA[data_histo_wider$MMA==""] <- 0
data_histo_wider$MMA <- as.numeric(data_histo_wider$MMA)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Hydropic vacuolation'] <- 'hydropic_vacuolation'
data_histo_wider$hydropic_vacuolation[data_histo_wider$hydropic_vacuolation == "Hydropic vacuolation"] <- 1
data_histo_wider$hydropic_vacuolation[data_histo_wider$hydropic_vacuolation==""] <- 0
data_histo_wider$hydropic_vacuolation <- as.numeric(data_histo_wider$hydropic_vacuolation)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Lipidosis'] <- 'lipidosis'
data_histo_wider$lipidosis[data_histo_wider$lipidosis == "Lipidosis"] <- 1
data_histo_wider$lipidosis[data_histo_wider$lipidosis==""] <- 0
data_histo_wider$lipidosis <- as.numeric(data_histo_wider$lipidosis)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Basophilic foci'] <- 'basophilic_foci'
data_histo_wider$basophilic_foci[data_histo_wider$basophilic_foci == "Basophilic foci"] <- 1
data_histo_wider$basophilic_foci[data_histo_wider$basophilic_foci==""] <- 0
data_histo_wider$basophilic_foci <- as.numeric(data_histo_wider$basophilic_foci)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Lymphocytic infiltration'] <- 'lymphocytic_infiltration'
data_histo_wider$lymphocytic_infiltration[data_histo_wider$lymphocytic_infiltration == "Lymphocytic infiltration"] <- 1
data_histo_wider$lymphocytic_infiltration[data_histo_wider$lymphocytic_infiltration==""] <- 0
data_histo_wider$lymphocytic_infiltration <- as.numeric(data_histo_wider$lymphocytic_infiltration)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Granuloma & Abscess'] <- 'granuloma_abscess'
data_histo_wider$granuloma_abscess[data_histo_wider$granuloma_abscess == "Granuloma & Abscess"] <- 1
data_histo_wider$granuloma_abscess[data_histo_wider$granuloma_abscess==""] <- 0
data_histo_wider$granuloma_abscess <- as.numeric(data_histo_wider$granuloma_abscess)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Clear cell foci'] <- 'clear_cell_foci'
data_histo_wider$clear_cell_foci[data_histo_wider$clear_cell_foci == "Clear cell foci"] <- 1
data_histo_wider$clear_cell_foci[data_histo_wider$clear_cell_foci==""] <- 0
data_histo_wider$clear_cell_foci <- as.numeric(data_histo_wider$clear_cell_foci)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Necrosis'] <- 'necrosis'
data_histo_wider$necrosis[data_histo_wider$necrosis == "Necrosis"] <- 1
data_histo_wider$necrosis[data_histo_wider$necrosis==""] <- 0
data_histo_wider$necrosis <- as.numeric(data_histo_wider$necrosis)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Eosinophilic foci'] <- 'eosinophilic_foci'
data_histo_wider$eosinophilic_foci[data_histo_wider$eosinophilic_foci == "Eosinophilic foci"] <- 1
data_histo_wider$eosinophilic_foci[data_histo_wider$eosinophilic_foci ==""] <- 0
data_histo_wider$eosinophilic_foci <- as.numeric(data_histo_wider$eosinophilic_foci)

colnames(data_histo_wider)[colnames(data_histo_wider) == 'Phospholipidosis'] <- 'phospholipidosis'
data_histo_wider$phospholipidosis[data_histo_wider$phospholipidosis == "Phospholipidosis"] <- 1
data_histo_wider$phospholipidosis[data_histo_wider$phospholipidosis ==""] <- 0
data_histo_wider$phospholipidosis <- as.numeric(data_histo_wider$phospholipidosis)

data_histo_lesions <- data_histo_wider %>%
  group_by(station, season, sex) %>%
  summarise(lymphocytic_infiltration = sum(lymphocytic_infiltration),
            MMA = sum(MMA),
            hydropic_vacuolation = sum(hydropic_vacuolation),
            necrosis = sum(necrosis),
            clear_cell_foci = sum(clear_cell_foci),
            lipidosis = sum(lipidosis),
            granuloma_abscess = sum(granuloma_abscess),
            basophilic_foci = sum(basophilic_foci),
            eosinophilic_foci = sum(eosinophilic_foci),
            phospholipidosis = sum(phospholipidosis))

data_histo_lesions <- merge(data_histo_lesions, table_station_season_sex, by.x=c('station', 'season', 'sex'), by.y=c('station', 'season', 'sex'))
data_histo_lesions <- data_histo_lesions %>% relocate(count_samples, .after = sex)

data_histo_lesions_prop <- data_histo_lesions %>%
  mutate(prop_lymphocytic_infiltration = lymphocytic_infiltration*100/count_samples,
         prop_MMA = MMA*100/count_samples,
         prop_hydropic_vacuolation = hydropic_vacuolation*100/count_samples,
         prop_necrosis = necrosis*100/count_samples,
         prop_clear_cell_foci = clear_cell_foci*100/count_samples,
         prop_lipidosis = lipidosis*100/count_samples,
         prop_granuloma_abscess = granuloma_abscess*100/count_samples,
         prop_basophilic_foci = basophilic_foci*100/count_samples,
         prop_eosinophilic_foci = eosinophilic_foci*100/count_samples,
         prop_phospholipidosis = phospholipidosis*100/count_samples) %>%
  select(-c(count_samples, lymphocytic_infiltration, MMA, hydropic_vacuolation, necrosis, clear_cell_foci, lipidosis, granuloma_abscess, basophilic_foci, eosinophilic_foci, phospholipidosis))

write_csv(data_histo_lesions_prop, "/home/valentin/Desktop/hepatovac/results/data_histo_lesions_prop.csv")

data_histo_lesions_prop$station <- as.factor(data_histo_lesions_prop$station) 
data_histo_lesions_prop$season <- as.factor(data_histo_lesions_prop$season) 
data_histo_lesions_prop$sex <- as.factor(data_histo_lesions_prop$sex) 

kruskal.test(prop_lymphocytic_infiltration ~ station, data = data_histo_lesions_prop) 
kruskal.test(prop_lymphocytic_infiltration ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_lymphocytic_infiltration ~ sex, data = data_histo_lesions_prop) 

kruskal.test(prop_MMA ~ station, data = data_histo_lesions_prop) # Significant difference detected in prevalence of MMA between stations, with higher prevalence in individuals from Seine compared to Canche.
kruskal.test(prop_MMA ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_MMA ~ sex, data = data_histo_lesions_prop) 
dunnTest(prop_MMA ~ station, data = data_histo_lesions_prop, method="bonferroni")
ggplot(data_histo_lesions_prop) +
  aes(x = season, y = prop_MMA, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_hydropic_vacuolation ~ station, data = data_histo_lesions_prop) # Significant difference detected in prevalence of hydropic vacuolation between stations, with higher prevalence in individuals from Seine compared to Canche.
kruskal.test(prop_hydropic_vacuolation ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_hydropic_vacuolation ~ sex, data = data_histo_lesions_prop) 
dunnTest(prop_hydropic_vacuolation ~ station, data = data_histo_lesions_prop, method="bonferroni")
ggplot(data_histo_lesions_prop) +
  aes(x = season, y = prop_hydropic_vacuolation, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_necrosis ~ station, data = data_histo_lesions_prop) # Significant difference detected in prevalence of necrosis between stations, with higher prevalence in individuals from Seine compared to Canche.
kruskal.test(prop_necrosis ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_necrosis ~ sex, data = data_histo_lesions_prop) 
dunnTest(prop_necrosis ~ station, data = data_histo_lesions_prop, method="bonferroni")
ggplot(data_histo_lesions_prop) +
  aes(x = season, y = prop_necrosis, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_clear_cell_foci ~ station, data = data_histo_lesions_prop)
kruskal.test(prop_clear_cell_foci ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_clear_cell_foci ~ sex, data = data_histo_lesions_prop) 

kruskal.test(prop_lipidosis ~ station, data = data_histo_lesions_prop)
kruskal.test(prop_lipidosis ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_lipidosis ~ sex, data = data_histo_lesions_prop) # Significant difference in prevalence of lipidosis between males and females, with males having higher prevalence of lipidosis compared to females.
dunnTest(prop_lipidosis ~ sex, data = data_histo_lesions_prop, method="bonferroni")
ggplot(data_histo_lesions_prop) +
  aes(x = season, y = prop_lipidosis, color = sex) +
  geom_point() +
  facet_wrap(~station)

kruskal.test(prop_granuloma_abscess ~ station, data = data_histo_lesions_prop) 
kruskal.test(prop_granuloma_abscess ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_granuloma_abscess ~ sex, data = data_histo_lesions_prop) 

kruskal.test(prop_basophilic_foci ~ station, data = data_histo_lesions_prop) 
kruskal.test(prop_basophilic_foci ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_basophilic_foci ~ sex, data = data_histo_lesions_prop) 

kruskal.test(prop_eosinophilic_foci ~ station, data = data_histo_lesions_prop) 
kruskal.test(prop_eosinophilic_foci ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_eosinophilic_foci ~ sex, data = data_histo_lesions_prop) 

kruskal.test(prop_phospholipidosis ~ station, data = data_histo_lesions_prop) 
kruskal.test(prop_phospholipidosis ~ season, data = data_histo_lesions_prop)
kruskal.test(prop_phospholipidosis ~ sex, data = data_histo_lesions_prop) 

# Test if there is a relationship between the prevalence of each lesions and log_prop_hepatovac_area 
priors <- get_prior(log_prop_hepatovac_area ~ 
                      lymphocytic_infiltration + 
                      MMA +
                      hydropic_vacuolation +
                      necrosis +
                      clear_cell_foci +
                      lipidosis +
                      granuloma_abscess +
                      basophilic_foci +
                      eosinophilic_foci +
                      phospholipidosis +
                      #age_month + 
                      sex +
                      season + 
                      station,
                    data = data_histo_wider, family = gaussian())

brm_histo_lesions_model <- brms::brm(log_prop_hepatovac_area ~ 
                                       lymphocytic_infiltration + 
                                       MMA +
                                       hydropic_vacuolation +
                                       necrosis +
                                       clear_cell_foci +
                                       lipidosis +
                                       granuloma_abscess +
                                       basophilic_foci +
                                       eosinophilic_foci +
                                       phospholipidosis +
                                       #age_month + 
                                       sex +
                                       season + 
                                       station, 
                                     data=data_histo_wider, 
                                     prior = priors)

summary(brm_histo_lesions_model)

## plot_histopathology -------------------------------------------------
(plot_histopathology <- ggarrange(plot_categories, plot_lesions, ncol = 2, nrow = 1, labels = c("A","B")))

ggsave("plot_histopathology.png", 
       plot = plot_histopathology, 
       width = 800, 
       height = 600, 
       dpi = 300, 
       units = "mm", 
       bg = "white")

rm(list = ls()[! ls() %in% c("data_TOXEM")])

## biomarker_data_TOXEM ----------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### AchE_µmol_min_mg_prot ***--------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(AchE_µmol_min_mg_prot) & !is.na(age_month))

ggplot(data_TOXEM) +
  aes(x = season, y = AchE_µmol_min_mg_prot, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_TOXEM %>%
  group_by(station, season) %>%
  summarise(mean_AchE_µmol_min_mg_prot = mean(AchE_µmol_min_mg_prot),
            sd_AchE_µmol_min_mg_prot = sd(AchE_µmol_min_mg_prot))

aov01 <- aov(AchE_µmol_min_mg_prot ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in AchE_µmol_min_mg_prot level between fishes from Seine and Canche, between fishes sampled in Winter vs Summer but not between females and males.

# Linear regression
ggplot(data_TOXEM) +
  aes(AchE_µmol_min_mg_prot, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = F) #+ facet_wrap(~ station)

cor.test(x = data_TOXEM$AchE_µmol_min_mg_prot, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# statistically significant positive correlation

lm_AchE <- lm(log_prop_hepatovac_area ~ AchE_µmol_min_mg_prot, data = data_TOXEM)
summary(lm_AchE)
# lm_AchE model show a significant effect on log_prop_hepatovac_area.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(AchE_µmol_min_mg_prot, log_log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_AchE <- gam(log_prop_hepatovac_area ~ s(AchE_µmol_min_mg_prot) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_AchE)
# gam_AchE model shows a significant effect on log_prop_hepatovac_area.

# Polynomial
poly2_AchE <- lm(log_prop_hepatovac_area ~ poly(AchE_µmol_min_mg_prot, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_AchE)
# poly2_AchE model show a significant effect on the response

poly3_AchE <- lm(log_prop_hepatovac_area ~ poly(AchE_µmol_min_mg_prot, 3) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(poly3_AchE)
# poly3_AchE model show a significant effect on the response

(plot_AchE <- ggplot(data_TOXEM) +
  aes(AchE_µmol_min_mg_prot, log_prop_hepatovac_area) + #, col = season, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 1), se = F)) #+ facet_wrap(~ station)

# Compare models
vif(lm_AchE)
vif(poly2_AchE)
vif(poly3_AchE)
# vif < 10 -> no collinearity between variables

summary(lm_AchE)$r.squared
summary(gam_AchE)$r.sq
summary(poly2_AchE)$r.squared
summary(poly3_AchE)$r.squared

AIC(lm_AchE, poly2_AchE)
BIC(lm_AchE, poly2_AchE)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_AchE, poly2_AchE, test = "F")

# Plot fitted values
plot(data_TOXEM$AchE_µmol_min_mg_prot, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
#points(data_TOXEM$AchE_µmol_min_mg_prot, fitted(lm_AchE), col = "green", pch = 4)
points(data_TOXEM$AchE_µmol_min_mg_prot, fitted(lm_AchE), col = "blue", pch = 3)
points(data_TOXEM$AchE_µmol_min_mg_prot, fitted(poly2_AchE), col = "red", pch = 2)
#legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$AchE_µmol_min_mg_prot, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$AchE_µmol_min_mg_prot, predict(lm_AchE), col = "green")
points(data_TOXEM$AchE_µmol_min_mg_prot, predict(gam_AchE), col = "blue")
points(data_TOXEM$AchE_µmol_min_mg_prot, predict(poly2_AchE), col = "red")

# Check assumptions
plot(lm_AchE)
res <- resid(poly2_AchE)
plot(fitted(poly2_AchE), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly(AchE_µmol_min_mg_prot, 2) model show a significant effect of AchE_µmol_min_mg_prot on log_prop_hepatovac_area.

### comet_percentage_DNA_tail ***--------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(comet_percentage_DNA_tail) & !is.na(age_month))

ggplot(data_TOXEM) +
  aes(x = season, y = comet_percentage_DNA_tail, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

ggplot(data_TOXEM) +
  aes(length_cm, comet_percentage_DNA_tail) +
  geom_point() +
  geom_smooth(method = lm, se = F)

cor.test(y = data_clean$comet_percentage_DNA_tail, x = data_clean$length_cm, method = "pearson", use = "complete.obs")

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_comet_percentage_DNA_tail = mean(comet_percentage_DNA_tail),
            sd_comet_percentage_DNA_tail = sd(comet_percentage_DNA_tail))

aov01 <- aov(comet_percentage_DNA_tail ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in comet_percentage_CNA_tail level between fishes from Seine and Canche, between fishes sampled in Winter vs Summer but not between females and males.

# Linear regression
(plot_comet_percentage_DNA_tail <- ggplot(data_TOXEM) +
  aes(comet_percentage_DNA_tail, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station))

cor.test(x = data_clean$comet_percentage_DNA_tail, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# statistically significant negative correlation

lm_comet <- lm(log_prop_hepatovac_area ~ comet_percentage_DNA_tail + poly(age_month, 2), data = data_clean)
summary(lm_comet)
# comet_percentage_DNA_tail model shows a significant effect on the response variable.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(comet_percentage_DNA_tail, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_comet <- gam(log_prop_hepatovac_area ~ s(comet_percentage_DNA_tail) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_comet)
# comet_percentage_DNA_tail model shows a significant effect on the response variable.

# Polynomial
poly2_comet <- lm(log_prop_hepatovac_area ~ poly(comet_percentage_DNA_tail, 2) + season + poly(age_month, 2), data = data_clean)
summary(poly2_comet)
# comet_percentage_DNA_tail model shows a significant effect on the response variable.

poly3_comet <- lm(log_prop_hepatovac_area ~ poly(comet_percentage_DNA_tail, 3) + season + poly(age_month, 2), data = data_clean)
summary(poly3_comet)
# comet_percentage_DNA_tail model shows a significant effect on the response variable.

ggplot(data_clean) +
  aes(comet_percentage_DNA_tail, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)


# Compare models
vif(lm_comet)
vif(gam_comet)
vif(poly2_comet)
vif(poly3_comet)
# vif < 10 -> no collinearity between variables

summary(lm_comet)$r.squared
summary(gam_comet)$r.sq
summary(poly2_comet)$r.squared
summary(poly3_comet)$r.squared

AIC(lm_comet, gam_comet, poly2_comet, poly3_comet)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_comet, gam_comet, poly2_comet, poly3_comet, test = "F")

# Plot fitted values
plot(data_TOXEM$comet_percentage_DNA_tail, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$comet_percentage_DNA_tail, fitted(lm_comet), col = "green", pch = 4)
points(data_clean$comet_percentage_DNA_tail, fitted(gam_comet), col = "blue", pch = 3)
points(data_clean$comet_percentage_DNA_tail, fitted(poly2_comet), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$comet_percentage_DNA_tail, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$comet_percentage_DNA_tail, predict(lm_comet), col = "green")
points(data_TOXEM$comet_percentage_DNA_tail, predict(gam_comet), col = "blue")
points(data_TOXEM$comet_percentage_DNA_tail, predict(poly2_comet), col = "red")

# Check assumptions
plot(lm_comet)
res <- resid(lm_comet)
plot(fitted(lm_comet), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly2_comet model show a significant effect of comet_percentage_DNA_tail on log_prop_hepatovac_area.
# poly2_comet


### EROD_pmol_min_mg_prot --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(EROD_pmol_min_mg_prot))

ggplot(data_TOXEM) +
  aes(x = season, y = EROD_pmol_min_mg_prot, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(EROD_pmol_min_mg_prot ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in EROD_pmol_min_mg_prot level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(EROD_pmol_min_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_clean$EROD_pmol_min_mg_prot, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_EROD <- lm(log_prop_hepatovac_area ~ EROD_pmol_min_mg_prot + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_EROD)
# EROD_pmol_min_mg_prot does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(EROD_pmol_min_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_EROD <- gam(log_prop_hepatovac_area ~ s(EROD_pmol_min_mg_prot) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_EROD)
# EROD_pmol_min_mg_prot does not significantly affect the response in this model.

# Polynomial
poly2_EROD <- lm(log_prop_hepatovac_area ~ poly(EROD_pmol_min_mg_prot, 2) + season + poly(age_month, 2), data = data_clean)
summary(poly2_EROD)
# EROD_pmol_min_mg_prot significantly affect the response in this model.

poly3_EROD <- lm(log_prop_hepatovac_area ~ poly(EROD_pmol_min_mg_prot, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_EROD)
# EROD_pmol_min_mg_prot does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(EROD_pmol_min_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

### TBARS_nmol_eq_MDA_mg_prot --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(TBARS_nmol_eq_MDA_mg_prot))

ggplot(data_TOXEM) +
  aes(x = season, y = TBARS_nmol_eq_MDA_mg_prot, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_TBARS_nmol_eq_MDA_mg_prot = mean(TBARS_nmol_eq_MDA_mg_prot), 
            sd_TBARS_nmol_eq_MDA_mg_prot = sd(TBARS_nmol_eq_MDA_mg_prot))

aov01 <- aov(TBARS_nmol_eq_MDA_mg_prot ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in TBARS_nmol_eq_MDA_mg_prot level between fishes sampled in Winter vs Summer. No significant difference between females and males.

# Linear regression
ggplot(data_clean) +
  aes(TBARS_nmol_eq_MDA_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + facet_wrap(~ station)

cor.test(x = data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_TBARS <- lm(log_prop_hepatovac_area ~ TBARS_nmol_eq_MDA_mg_prot + season + poly(age_month, 2), data = data_TOXEM)
summary(lm_TBARS)
# TBARS_nmol_eq_MDA_mg_prot does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(TBARS_nmol_eq_MDA_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_TBARS <- gam(log_prop_hepatovac_area ~ s(TBARS_nmol_eq_MDA_mg_prot) + season + poly(age_month, 2), data = data_clean)
summary(gam_TBARS)
# TBARS_nmol_eq_MDA_mg_prot does not significantly affect the response in this model.

# Polynomial
poly2_TBARS <- lm(log_prop_hepatovac_area ~ poly(TBARS_nmol_eq_MDA_mg_prot, 2) + season + poly(age_month, 2), data = data_clean)
summary(poly2_TBARS)
# TBARS_nmol_eq_MDA_mg_prot does not significantly affect the response in this model.

poly3_TBARS <- lm(log_prop_hepatovac_area ~ poly(TBARS_nmol_eq_MDA_mg_prot, 3) + season + poly(age_month, 2), data = data_clean)
summary(poly3_TBARS)
# TBARS_nmol_eq_MDA_mg_prot does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(TBARS_nmol_eq_MDA_mg_prot, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

# Compare models
vif(lm_TBARS)
vif(gam_TBARS)
vif(poly2_TBARS)
vif(poly3_TBARS)
# vif < 10 -> no collinearity between variables

summary(lm_TBARS)$r.squared
summary(gam_TBARS)$r.sq
summary(poly2_TBARS)$r.squared
summary(poly3_TBARS)$r.squared

AIC(lm_TBARS, gam_TBARS, poly2_TBARS, poly3_TBARS)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_TBARS, gam_TBARS, poly2_TBARS, poly3_TBARS, test = "F")

# Plot fitted values
plot(data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$TBARS_nmol_eq_MDA_mg_prot, fitted(lm_TBARS), col = "green", pch = 4)
points(data_clean$TBARS_nmol_eq_MDA_mg_prot, fitted(gam_TBARS), col = "blue", pch = 3)
points(data_clean$TBARS_nmol_eq_MDA_mg_prot, fitted(poly2_TBARS), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, predict(lm_TBARS), col = "green")
points(data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, predict(gam_TBARS), col = "blue")
points(data_TOXEM$TBARS_nmol_eq_MDA_mg_prot, predict(poly2_TBARS), col = "red")

# Check assumptions
plot(lm_TBARS)
res <- resid(lm_TBARS)
plot(fitted(lm_TBARS), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### prot_carbo_nmol_mg --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(prot_carbo_nmol_mg))

ggplot(data_TOXEM) +
  aes(x = season, y = prot_carbo_nmol_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_prot_carbo_nmol_mg = mean(prot_carbo_nmol_mg), 
            sd_prot_carbo_nmol_mg = sd(prot_carbo_nmol_mg))

aov01 <- aov(prot_carbo_nmol_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in prot_carbo_nmol_mg level between fishes sampled in Winter vs Summer.

# Linear regression
ggplot(data_TOXEM) +
  aes(prot_carbo_nmol_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$prot_carbo_nmol_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_prot_carbo <- lm(log_prop_hepatovac_area ~ prot_carbo_nmol_mg + station * season * sex + poly(age_month, 2), data = data_clean)
summary(lm_prot_carbo)
# prot_carbo_nmol_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(prot_carbo_nmol_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_prot_carbo <- gam(log_prop_hepatovac_area ~ s(prot_carbo_nmol_mg) + station * season + poly(age_month, 2), data = data_clean)
summary(gam_prot_carbo)
# prot_carbo_nmol_mg does not significantly affect the response in this model.

# Polynomial
poly2_prot_carbo <- lm(log_prop_hepatovac_area ~ poly(prot_carbo_nmol_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_prot_carbo)
# prot_carbo_nmol_mg does not significantly affect the response in this model.

poly3_prot_carbo <- lm(log_prop_hepatovac_area ~ poly(prot_carbo_nmol_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_prot_carbo)
# prot_carbo_nmol_mg does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(prot_carbo_nmol_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)


### G6PDH_liver_IU_mg ***--------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(G6PDH_liver_IU_mg) & !is.na(age_month))

ggplot(data_TOXEM) +
  aes(x = season, y = G6PDH_liver_IU_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_G6PDH_liver_IU_mg = mean(G6PDH_liver_IU_mg),
            sd_G6PDH_liver_IU_mg = sd(G6PDH_liver_IU_mg))

aov01 <- aov(G6PDH_liver_IU_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in G6PDH_liver_IU_mg level between station and season. No significant difference between sex. Level of G6PDH_liver_IU_mg being significantly higher in fish sampled in the Seine estuary in summer comparatively to the other groups. 

# Linear regression
ggplot(data_clean) +
  aes(G6PDH_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_clean$G6PDH_liver_IU_mg, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_G6PDH <- lm(log_prop_hepatovac_area ~ G6PDH_liver_IU_mg + station * season + poly(age_month, 2), data = data_clean)
summary(lm_G6PDH)
# G6PDH_liver_IU_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(G6PDH_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_G6PDH <- gam(log_prop_hepatovac_area ~ s(G6PDH_liver_IU_mg) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_G6PDH)
# G6PDH_liver_IU_mg does not significantly affect the response in this model.

# Polynomial
poly2_G6PDH <- lm(log_prop_hepatovac_area ~ poly(G6PDH_liver_IU_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_G6PDH)
# poly2_G6PDH model shows a significant effect on the response.

poly3_G6PDH <- lm(log_prop_hepatovac_area ~ poly(G6PDH_liver_IU_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_G6PDH)
# poly3_G6PDH model shows a significant effect on the response.

(plot_G6PDH <- ggplot(data_clean) +
  aes(G6PDH_liver_IU_mg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F))

# Compare models
vif(lm_G6PDH)
vif(gam_G6PDH)
vif(poly2_G6PDH)
vif(poly3_G6PDH)
# vif < 10 -> no collinearity between variables

summary(lm_G6PDH)$r.squared
summary(gam_G6PDH)$r.sq
summary(poly2_G6PDH)$r.squared
summary(poly3_G6PDH)$r.squared

AIC(lm_G6PDH, gam_G6PDH, poly2_G6PDH, poly3_G6PDH)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_G6PDH, gam_G6PDH, poly2_G6PDH, poly3_G6PDH, test = "F")

# Plot fitted values
plot(data_TOXEM$G6PDH_liver_IU_mg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$G6PDH_liver_IU_mg, fitted(lm_G6PDH), col = "green", pch = 4)
points(data_clean$G6PDH_liver_IU_mg, fitted(gam_G6PDH), col = "blue", pch = 3)
points(data_clean$G6PDH_liver_IU_mg, fitted(poly2_G6PDH), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$G6PDH_liver_IU_mg, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$G6PDH_liver_IU_mg, predict(lm_G6PDH), col = "green")
points(data_TOXEM$G6PDH_liver_IU_mg, predict(gam_G6PDH), col = "blue")
points(data_TOXEM$G6PDH_liver_IU_mg, predict(poly2_G6PDH), col = "red")

# Check assumptions
plot(lm_G6PDH)
res <- resid(lm_G6PDH)
plot(fitted(lm_G6PDH), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly3_G6PDH model show a significant effect of G6PDH_liver_IU_mg on log_prop_hepatovac_area.
# poly3_G6PDH

### CS_liver_IU_mg --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(CS_liver_IU_mg))

ggplot(data_TOXEM) +
  aes(x = season, y = CS_liver_IU_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(CS_liver_IU_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in citrate synthase (IU.mg^-1) level between station and season. No significant difference between sex. 
# Fish from the Canche estuary present a significantly higher level of CS_liver_IU_mg than fish from the Seine estuary.
# Level of G6PDH_liver_IU_mg being significantly higher in fish sampled in the Seine estuary in summer comparatively to the other groups. Fish from the Canche estuary sampled in summer shows a significantly higher level of CS_liver_IU_mg than fish from the Canche estuary sampled in winter and fish from the Seine estuary. 

# Linear regression
ggplot(data_TOXEM) +
  aes(CS_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_clean$CS_liver_IU_mg, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_CS_liver <- lm(log_prop_hepatovac_area ~ CS_liver_IU_mg + station * season * sex + poly(age_month, 2), data = data_clean)
summary(lm_CS_liver)
# CS_liver_IU_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(CS_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_CS_liver <- gam(log_prop_hepatovac_area ~ s(CS_liver_IU_mg) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_CS_liver)
# CS_liver_IU_mg does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(CS_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_CS_liver <- lm(log_prop_hepatovac_area ~ poly(CS_liver_IU_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_CS_liver)
# CS_liver_IU_mg does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(CS_liver_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_CS_liver <- lm(log_prop_hepatovac_area ~ poly(CS_liver_IU_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_CS_liver)
# CS_liver_IU_mg does not significantly affect the response in this model.

### CS_muscle_IU_mg --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(CS_muscle_IU_mg))
#data_clean <- na.omit(data_TOXEM[, c("CS_muscle_IU_mg", "log_prop_hepatovac_area", "HSI", "station", "season","sex")])

ggplot(data_TOXEM) +
  aes(x = season, y = CS_muscle_IU_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

aov01 <- aov(CS_muscle_IU_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

# Linear regression
ggplot(data_clean) +
  aes(CS_muscle_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_clean$CS_muscle_IU_mg, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_CS_muscle <- lm(log_prop_hepatovac_area ~ CS_muscle_IU_mg + station * season * sex + poly(age_month, 2), data = data_clean)
summary(lm_CS_muscle)
# CS_muscle_IU_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(CS_muscle_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_CS_muscle <- gam(log_prop_hepatovac_area ~ s(CS_muscle_IU_mg) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_CS_muscle)
# CS_muscle_IU_mg does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(CS_muscle_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_CS_muscle <- lm(log_prop_hepatovac_area ~ poly(CS_muscle_IU_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_CS_muscle)
# CS_muscle_IU_mg does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(CS_muscle_IU_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_CS_muscle <- lm(log_prop_hepatovac_area ~ poly(CS_muscle_IU_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_CS_muscle)
# CS_muscle_IU_mg does not significantly affect the response in this model.

# Compare models
vif(lm_CS_muscle)
vif(gam_CS_muscle)
vif(poly2_CS_muscle)
vif(poly3_CS_muscle)
# vif < 10 -> no collinearity between variables

summary(lm_CS_muscle)$r.squared
summary(gam_CS_muscle)$r.sq
summary(poly2_CS_muscle)$r.squared
summary(poly3_CS_muscle)$r.squared

AIC(lm_CS_muscle, gam_CS_muscle, poly2_CS_muscle, poly3_CS_muscle)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_CS_muscle, gam_CS_muscle, poly2_CS_muscle, poly3_CS_muscle, test = "F")

# Plot fitted values
plot(data_TOXEM$CS_muscle_IU_mg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$CS_muscle_IU_mg, fitted(lm_CS_muscle), col = "green", pch = 4)
points(data_clean$CS_muscle_IU_mg, fitted(gam_CS_muscle), col = "blue", pch = 3)
points(data_clean$CS_muscle_IU_mg, fitted(poly2_CS_muscle), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$CS_muscle_IU_mg, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$CS_muscle_IU_mg, predict(lm_CS_muscle), col = "green")
points(data_TOXEM$CS_muscle_IU_mg, predict(gam_CS_muscle), col = "blue")
points(data_TOXEM$CS_muscle_IU_mg, predict(poly2_CS_muscle), col = "red")

# Check assumptions
plot(lm_CS_muscle)
res <- resid(lm_CS_muscle)
plot(fitted(lm_CS_muscle), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### prot_liver_mg_prot_mL_enzymes --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(prot_liver_mg_prot_mL_enzymes))
#data_clean <- na.omit(data_TOXEM[, c("prot_liver_mg_prot_mL_enzymes", "log_prop_hepatovac_area", "HSI", "station", "season","sex")])

aov01 <- aov(prot_liver_mg_prot_mL_enzymes ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

# Linear regression
ggplot(data_clean) +
  aes(prot_liver_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_clean$prot_liver_mg_prot_mL_enzymes, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_prot_liver <- lm(log_prop_hepatovac_area ~ prot_liver_mg_prot_mL_enzymes + station * season * sex + poly(age_month, 2), data = data_clean)
summary(lm_prot_liver)
# prot_liver_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(prot_liver_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_prot_liver <- gam(log_prop_hepatovac_area ~ s(prot_liver_mg_prot_mL_enzymes) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_prot_liver)
# prot_liver_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(prot_liver_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_prot_liver <- lm(log_prop_hepatovac_area ~ poly(prot_liver_mg_prot_mL_enzymes, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_prot_liver)
# prot_liver_mg_prot_mL_enzymes does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(prot_liver_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_prot_liver <- lm(log_prop_hepatovac_area ~ poly(prot_liver_mg_prot_mL_enzymes, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_prot_liver)
# prot_liver_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Compare models
vif(lm_prot_liver)
vif(gam_prot_liver)
vif(poly2_prot_liver)
vif(poly3_prot_liver)
# vif < 10 -> no collinearity between variables

summary(lm_prot_liver)$r.squared
summary(gam_prot_liver)$r.sq
summary(poly2_prot_liver)$r.squared
summary(poly3_prot_liver)$r.squared

AIC(lm_prot_liver, gam_prot_liver, poly2_prot_liver, poly3_prot_liver)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_prot_liver, gam_prot_liver, poly2_prot_liver, poly3_prot_liver, test = "F")

# Plot fitted values
plot(data_TOXEM$prot_liver_mg_prot_mL_enzymes, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$prot_liver_mg_prot_mL_enzymes, fitted(lm_prot_liver), col = "green", pch = 4)
points(data_clean$prot_liver_mg_prot_mL_enzymes, fitted(gam_prot_liver), col = "blue", pch = 3)
points(data_clean$prot_liver_mg_prot_mL_enzymes, fitted(poly2_prot_liver), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$prot_liver_mg_prot_mL_enzymes, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$prot_liver_mg_prot_mL_enzymes, predict(lm_prot_liver), col = "green")
points(data_TOXEM$prot_liver_mg_prot_mL_enzymes, predict(gam_prot_liver), col = "blue")
points(data_TOXEM$prot_liver_mg_prot_mL_enzymes, predict(poly2_prot_liver), col = "red")

# Check assumptions
plot(lm_prot_liver)
res <- resid(lm_prot_liver)
plot(fitted(lm_prot_liver), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

### prot_muscle_mg_prot_mL_enzymes --------------------------------------------------------------
# Remove rows with missing values in relevant columns
data_clean <- subset(data_TOXEM, !is.na(prot_muscle_mg_prot_mL_enzymes))
#data_clean <- na.omit(data_TOXEM[, c("prot_muscle_mg_prot_mL_enzymes", "log_prop_hepatovac_area", "HSI", "station", "season","sex")])

aov01 <- aov(prot_muscle_mg_prot_mL_enzymes ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

# Linear regression
ggplot(data_clean) +
  aes(prot_muscle_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_clean$prot_muscle_mg_prot_mL_enzymes, y = data_clean$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_prot_muscle <- lm(log_prop_hepatovac_area ~ prot_muscle_mg_prot_mL_enzymes + station * season * sex + poly(age_month, 2), data = data_clean)
summary(lm_prot_muscle)
# prot_muscle_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_clean) +
  aes(prot_muscle_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_prot_muscle <- gam(log_prop_hepatovac_area ~ s(prot_muscle_mg_prot_mL_enzymes) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_prot_muscle)
# prot_muscle_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Polynomial
ggplot(data_clean) +
  aes(prot_muscle_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F) + facet_wrap(~ station)

poly2_prot_muscle <- lm(log_prop_hepatovac_area ~ poly(prot_muscle_mg_prot_mL_enzymes, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_prot_muscle)
# prot_muscle_mg_prot_mL_enzymes does not significantly affect the response in this model.

ggplot(data_clean) +
  aes(prot_muscle_mg_prot_mL_enzymes, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F) + facet_wrap(~ station)

poly3_prot_muscle <- lm(log_prop_hepatovac_area ~ poly(prot_muscle_mg_prot_mL_enzymes, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_prot_muscle)
# prot_muscle_mg_prot_mL_enzymes does not significantly affect the response in this model.

# Compare models
vif(lm_prot_muscle)
vif(gam_prot_muscle)
vif(poly2_prot_muscle)
vif(poly3_prot_muscle)
# vif < 10 -> no collinearity between variables

summary(lm_prot_muscle)$r.squared
summary(gam_prot_muscle)$r.sq
summary(poly2_prot_muscle)$r.squared
summary(poly3_prot_muscle)$r.squared

AIC(lm_prot_muscle, gam_prot_muscle, poly2_prot_muscle, poly3_prot_muscle)
# Lower values indicate better model fit, penalizing complexity.

anova(lm_prot_muscle, gam_prot_muscle, poly2_prot_muscle, poly3_prot_muscle, test = "F")

# Plot fitted values
plot(data_TOXEM$prot_muscle_mg_prot_mL_enzymes, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$prot_muscle_mg_prot_mL_enzymes, fitted(lm_prot_muscle), col = "green", pch = 4)
points(data_clean$prot_muscle_mg_prot_mL_enzymes, fitted(gam_prot_muscle), col = "blue", pch = 3)
points(data_clean$prot_muscle_mg_prot_mL_enzymes, fitted(poly2_prot_muscle), col = "red", pch = 2)
legend("topleft", legend = c("LM", "GAM", "PM"), col = c("green", "blue", "red"), pch = c(2, 3))

plot(data_TOXEM$prot_muscle_mg_prot_mL_enzymes, data_TOXEM$log_prop_hepatovac_area)
points(data_TOXEM$prot_muscle_mg_prot_mL_enzymes, predict(lm_prot_muscle), col = "green")
points(data_TOXEM$prot_muscle_mg_prot_mL_enzymes, predict(gam_prot_muscle), col = "blue")
points(data_TOXEM$prot_muscle_mg_prot_mL_enzymes, predict(poly2_prot_muscle), col = "red")

# Check assumptions
plot(lm_prot_muscle)
res <- resid(lm_prot_muscle)
plot(fitted(lm_prot_muscle), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

## plot_biomarker_data_TOXEM -------------------------------------------------
ggarrange(plot_AchE, plot_G6PDH, plot_comet_percentage_DNA_tail, ncol = 2, nrow = 2, labels = c("A","B","C"))

## lipids_data_TOXEM ------------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

## neutral_lipids ----------------------------------------------------------
### lipids_SE_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_SE_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_SE_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_SE_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_SE_µg_mg = mean(lipids_SE_µg_mg))

aov01 <- aov(lipids_SE_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_SE_µg_mg levels between fishes from the different groups.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_SE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_SE_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_SE <- lm(log_prop_hepatovac_area ~ lipids_SE_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_SE)
# lipids_SE_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_SE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_SE <- gam(log_prop_hepatovac_area ~ s(lipids_SE_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_SE)
# lipids_SE_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_SE <- lm(log_prop_hepatovac_area ~ poly(lipids_SE_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_SE)
# lipids_SE_µg_mg does not significantly affect the response in this model.

poly3_lipids_SE <- lm(log_prop_hepatovac_area ~ poly(lipids_SE_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_SE)
# lipids_SE_µg_mg significantly affect the response in this model.


### lipids_GE_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_GE_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_GE_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_GE_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_lipids_GE_µg_mg = mean(lipids_GE_µg_mg),
            sd_lipids_GE_µg_mg = sd(lipids_GE_µg_mg))

aov01 <- aov(lipids_GE_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in lipids_GE_µg_mg levels between fishes sampled at different season.
# Fish sampled in summer present a significantly higher level of lipids_FS_µg_mg than fish sampled in winter. 
# No significant difference in lipids_FS_µg_mg levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_GE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_GE_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_GE <- lm(log_prop_hepatovac_area ~ lipids_GE_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_GE)
# lipids_GE_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_GE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_GE <- gam(log_prop_hepatovac_area ~ s(lipids_GE_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_GE)
# lipids_GE_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_GE <- lm(log_prop_hepatovac_area ~ poly(lipids_GE_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_GE)
# lipids_GE_µg_mg does not significantly affect the response in this model.

poly3_lipids_GE <- lm(log_prop_hepatovac_area ~ poly(lipids_GE_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_GE)
# lipids_GE_µg_mg does not significantly affect the response in this model.

### lipids_TG_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_TG_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_TG_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_TG_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_lipids_TG_µg_mg = mean(lipids_TG_µg_mg), 
            sd_lipids_TG_µg_mg = sd(lipids_TG_µg_mg))

aov01 <- aov(lipids_TG_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in lipids_TG_µg_mg levels between fishes sampled at the different station and season.
# Fish sampled in Seine present a significantly higher level of lipids_TG_µg_mg than fish sampled in Canche.
# Fish sampled in summer present a significantly higher level of lipids_TG_µg_mg than fish sampled in winter.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_TG_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_TG_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_TG <- lm(log_prop_hepatovac_area ~ lipids_TG_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_TG)
# lipids_TG_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_TG_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_TG <- gam(log_prop_hepatovac_area ~ s(lipids_TG_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_TG)
# lipids_TG_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_TG <- lm(log_prop_hepatovac_area ~ poly(lipids_TG_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_TG)
# lipids_TG_µg_mg does not significantly affect the response in this model.

poly3_lipids_TG <- lm(log_prop_hepatovac_area ~ poly(lipids_TG_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_TG)
# lipids_TG_µg_mg does not significantly affect the response in this model.

### lipids_FFA_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_FFA_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_FFA_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_FFA_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_FFA_µg_mg = mean(lipids_FFA_µg_mg))

aov01 <- aov(lipids_FFA_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_FFA_µg_mg levels between fishes from the different groups.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_FFA_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_FFA_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_FFA <- lm(log_prop_hepatovac_area ~ lipids_FFA_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_FFA)
# lipids_FFA_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_FFA_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_FFA <- gam(log_prop_hepatovac_area ~ s(lipids_FFA_µg_mg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_FFA)
# lipids_FFA_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_FFA <- lm(log_prop_hepatovac_area ~ poly(lipids_FFA_µg_mg, 2) + season + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_FFA)

poly3_lipids_FFA <- lm(log_prop_hepatovac_area ~ poly(lipids_FFA_µg_mg, 3) + season + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_FFA)

### lipids_ALC_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_ALC_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_ALC_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_ALC_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_ALC_µg_mg = mean(lipids_ALC_µg_mg))

aov01 <- aov(lipids_ALC_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_ALC_µg_mg levels between fishes from the different groups.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_ALC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_ALC_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_ALC <- lm(log_prop_hepatovac_area ~ lipids_ALC_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_ALC)
# lipids_ALC_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_ALC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_ALC <- gam(log_prop_hepatovac_area ~ s(lipids_ALC_µg_mg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_ALC)
# lipids_ALC_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_ALC <- lm(log_prop_hepatovac_area ~ poly(lipids_ALC_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_ALC)
# lipids_ALC_µg_mg does not significantly affect the response in this model.

poly3_lipids_ALC <- lm(log_prop_hepatovac_area ~ poly(lipids_ALC_µg_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_ALC)
# lipids_ALC_µg_mg does not significantly affect the response in this model.

### lipids_FS_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_FS_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_FS_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_FS_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_lipids_FS_µg_mg = mean(lipids_FS_µg_mg),
            sd_lipids_FS_µg_mg = sd(lipids_FS_µg_mg))

aov01 <- aov(lipids_FS_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in lipids_FS_µg_mg levels between fishes sampled at different season.
# Fish sampled in summer present a significantly higher level of lipids_FS_µg_mg than fish sampled in winter. 
# No significant difference in lipids_FS_µg_mg levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_FS_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_FS_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_FS <- lm(log_prop_hepatovac_area ~ lipids_FS_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_FS)
# lipids_FS_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_FS_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_FS <- gam(log_prop_hepatovac_area ~ s(lipids_FS_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_FS)
# lipids_FS_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_FS <- lm(log_prop_hepatovac_area ~ poly(lipids_FS_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_FS)
# lipids_FS_µg_mg does not significantly affect the response in this model.

poly3_lipids_FS <- lm(log_prop_hepatovac_area ~ poly(lipids_FS_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_FS)
# lipids_FS_µg_mg does not significantly affect the response in this model.

### total_neutral_lipids_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_neutral_lipids_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("total_neutral_lipids_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_neutral_lipids_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_total_neutral_lipids_µg_mg = mean(total_neutral_lipids_µg_mg),
            sd_total_neutral_lipids_µg_mg = sd(total_neutral_lipids_µg_mg))

aov01 <- aov(total_neutral_lipids_µg_mg ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_neutral_lipids_µg_mg levels between fishes from the different season.
# Fish sampled in winter present a significantly lower level of total_neutral_lipids_µg_mg than fish sampled in summer.
# No significant difference in total_neutral_lipids_µg_mg levels between sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_neutral_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_neutral_lipids_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_total_neutral <- lm(log_prop_hepatovac_area ~ total_neutral_lipids_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_neutral)
# total_neutral_lipids_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_neutral_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_total_neutral <- gam(log_prop_hepatovac_area ~ s(total_neutral_lipids_µg_mg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_neutral)
# total_neutral_lipids_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_total_neutral <- lm(log_prop_hepatovac_area ~ poly(total_neutral_lipids_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_total_neutral)
# total_neutral_lipids_µg_mg does not significantly affect the response in this model.

poly3_total_neutral <- lm(log_prop_hepatovac_area ~ poly(total_neutral_lipids_µg_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_total_neutral)
# total_neutral_lipids_µg_mg does not significantly affect the response in this model.

### total_reserve_lipids_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_reserve_lipids_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("total_reserve_lipids_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_reserve_lipids_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_total_reserve_lipids_µg_mg = mean(total_reserve_lipids_µg_mg), 
            sd_total_reserve_lipids_µg_mg = sd(total_reserve_lipids_µg_mg))

aov01 <- aov(total_reserve_lipids_µg_mg ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_reserve_lipids_µg_mg levels between season


# Linear regression
ggplot(data_TOXEM) +
  aes(total_reserve_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_reserve_lipids_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_total_reserve <- lm(log_prop_hepatovac_area ~ total_reserve_lipids_µg_mg + season + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_reserve)
# total_reserve_lipids_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_reserve_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_total_reserve <- gam(log_prop_hepatovac_area ~ s(total_reserve_lipids_µg_mg) + season + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_reserve)
# total_reserve_lipids_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_total_reserve <- lm(log_prop_hepatovac_area ~ poly(total_reserve_lipids_µg_mg, 2) + station + season + poly(age_month, 2), data = data_clean)
summary(poly2_total_reserve)
# total_reserve_lipids_µg_mg does not significantly affect the response in this model.

poly3_total_reserve <- lm(log_prop_hepatovac_area ~ poly(total_reserve_lipids_µg_mg, 3) + station + season + poly(age_month, 2), data = data_clean)
summary(poly3_total_reserve)
# total_reserve_lipids_µg_mg does not significantly affect the response in this model.

## polar_lipids ------------------------------------------------------------
### lipids_SPG_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_SPG_µg_mg) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("lipids_SPG_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_SPG_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_lipids_SPG_µg_mg = mean(lipids_SPG_µg_mg),
            sd_lipids_SPG_µg_mg = sd(lipids_SPG_µg_mg))

aov01 <- aov(lipids_SPG_µg_mg ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_SPG_µg_mg levels between fishes from the different groups.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_SPG_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_SPG_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_SPG <- lm(log_prop_hepatovac_area ~ lipids_SPG_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_SPG)
# lipids_SPG_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_SPG_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_SPG <- gam(log_prop_hepatovac_area ~ s(lipids_SPG_µg_mg) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(gam_lipids_SPG)
# lipids_SPG_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_SPG <- lm(log_prop_hepatovac_area ~ poly(lipids_SPG_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_SPG)
# lipids_SPG_µg_mg does not significantly affect the response in this model.

poly3_lipids_SPG <- lm(log_prop_hepatovac_area ~ poly(lipids_SPG_µg_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_SPG)
# lipids_SPG_µg_mg does not significantly affect the response in this model.

### lipids_LPC_µg_mg ***--------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_LPC_µg_mg) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("lipids_LPC_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_LPC_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_LPC_µg_mg = mean(lipids_LPC_µg_mg),
            sd_lipids_LPC_µg_mg = sd(lipids_LPC_µg_mg))

aov01 <- aov(lipids_LPC_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_LPC_µg_mg levels between fishes from the different groups 

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_LPC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_LPC_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_LPC <- lm(log_prop_hepatovac_area ~ lipids_LPC_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_LPC)
# lipids_LPC_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_LPC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_LPC <- gam(log_prop_hepatovac_area ~ s(lipids_LPC_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_LPC)
# lipids_LPC_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_LPC <- lm(log_prop_hepatovac_area ~ poly(lipids_LPC_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_LPC)
# lipids_LPC_µg_mg does not significantly affect the response in this model.

poly3_lipids_LPC <- lm(log_prop_hepatovac_area ~ poly(lipids_LPC_µg_mg, 3) + season, data = data_clean)
summary(poly3_lipids_LPC)
# lipids_LPC_µg_mg significantly affect the response in this model.

(plot_lipids_LPC <- ggplot(data_clean) +
  aes(lipids_LPC_µg_mg, log_prop_hepatovac_area) + #, col = station, shape = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F))

# Compare models
vif(poly3_lipids_LPC)
# vif < 10 -> no collinearity between variables

summary(poly3_lipids_LPC)$r.squared

AIC(poly3_lipids_LPC)
# Lower values indicate better model fit, penalizing complexity.

# Plot fitted values
plot(data_TOXEM$lipids_LPC_µg_mg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$lipids_LPC_µg_mg, fitted(poly3_lipids_LPC), col = "red", pch = 2)

# Check assumptions
plot(poly3_lipids_LPC)
res <- resid(poly3_lipids_LPC)
plot(fitted(poly3_lipids_LPC), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly3_lipids_LPC

### lipids_PC_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_PC_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_PC_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_PC_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_PC_µg_mg = mean(lipids_PC_µg_mg))

aov01 <- aov(lipids_PC_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_PC_µg_mg levels between fishes from the different groups 

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_PC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_PC_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_PC <- lm(log_prop_hepatovac_area ~ lipids_PC_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_PC)
# lipids_PC_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_PC_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_PC <- gam(log_prop_hepatovac_area ~ s(lipids_PC_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_PC)
# lipids_PC_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_PC <- lm(log_prop_hepatovac_area ~ poly(lipids_PC_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_PC)
# lipids_PC_µg_mg does not significantly affect the response in this model.

poly3_lipids_PC <- lm(log_prop_hepatovac_area ~ poly(lipids_PC_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_PC)
# lipids_PC_µg_mg does not significantly affect the response in this model.

### lipids_PS_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_PS_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_PS_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_PS_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_lipids_PS_µg_mg = mean(lipids_PS_µg_mg),
            sd_lipids_PS_µg_mg = sd(lipids_PS_µg_mg))

aov01 <- aov(lipids_PS_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in lipids_PS_µg_mg levels between fishes from the different season.
# Fish sampled in summer present a significantly higher level of lipids_PS_µg_mg than fish sampled in winter. 
# No significant difference in lipids_PE_µg_mg levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_PS_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_PS_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_PS <- lm(log_prop_hepatovac_area ~ lipids_PS_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_PS)
# lipids_PS_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_PS_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_PS <- gam(log_prop_hepatovac_area ~ s(lipids_PS_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_PS)
# lipids_PS_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_PS <- lm(log_prop_hepatovac_area ~ poly(lipids_PS_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_PS)
# lipids_PS_µg_mg does not significantly affect the response in this model.

poly3_lipids_PS <- lm(log_prop_hepatovac_area ~ poly(lipids_PS_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_PS)
# lipids_PS_µg_mg does not significantly affect the response in this model.

### lipids_PI_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_PI_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_PI_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_PI_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_PI_µg_mg = mean(lipids_PI_µg_mg))

aov01 <- aov(lipids_PI_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_PI_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_PI_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_PI <- lm(log_prop_hepatovac_area ~ lipids_PI_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_PI)
# lipids_PI_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_PI_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_PI <- gam(log_prop_hepatovac_area ~ s(lipids_PI_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_PI)
# lipids_PI_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_PI <- lm(log_prop_hepatovac_area ~ poly(lipids_PI_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_PI)
# lipids_PI_µg_mg does not significantly affect the response in this model.

poly3_lipids_PI <- lm(log_prop_hepatovac_area ~ poly(lipids_PI_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_PI)
# lipids_PI_µg_mg does not significantly affect the response in this model.

### lipids_CL_µg_mg ***--------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_CL_µg_mg) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("lipids_CL_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_CL_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_lipids_CL_µg_mg = mean(lipids_CL_µg_mg))

aov01 <- aov(lipids_CL_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in lipids_CL_µg_mg levels between fishes from the different groups 

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_CL_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_CL_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_CL <- lm(log_prop_hepatovac_area ~ lipids_CL_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_CL)
# lipids_CL_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_CL_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_CL <- gam(log_prop_hepatovac_area ~ s(lipids_CL_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_CL)
# lipids_CL_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_CL <- lm(log_prop_hepatovac_area ~ poly(lipids_CL_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_CL)
# lipids_CL_µg_mg does not significantly affect the response in this model.

poly3_lipids_CL <- lm(log_prop_hepatovac_area ~ poly(lipids_CL_µg_mg, 3) + season + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_CL)
# lipids_CL_µg_mg significantly affect the response in this model.

(plot_lipids_CL <- ggplot(data_clean) +
  aes(lipids_CL_µg_mg, log_prop_hepatovac_area) + #, col = station, shape = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F))

# Compare models
vif(poly3_lipids_CL)
# vif < 10 -> no collinearity between variables

summary(poly3_lipids_CL)$r.squared

# Plot fitted values
plot(data_TOXEM$lipids_CL_µg_mg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$lipids_CL_µg_mg, fitted(poly3_lipids_CL), col = "red", pch = 2)

# Check assumptions
plot(poly3_lipids_CL)
res <- resid(poly3_lipids_CL)
plot(fitted(poly3_lipids_CL), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# poly3_lipids_CL

### lipids_PE_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(lipids_PE_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("lipids_PE_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = lipids_PE_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_lipids_PE_µg_mg = mean(lipids_PE_µg_mg),
            sd_lipids_PE_µg_mg = sd(lipids_PE_µg_mg))

aov01 <- aov(lipids_PE_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in lipids_PE_µg_mg levels between fishes from the different season.
# Fish sampled in summer present a significantly higher level of lipids_PE_µg_mg than fish sampled in winter. 
# No significant difference in lipids_PE_µg_mg levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(lipids_PE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$lipids_PE_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_lipids_PE <- lm(log_prop_hepatovac_area ~ lipids_PE_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_lipids_PE)
# lipids_PE_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(lipids_PE_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_lipids_PE <- gam(log_prop_hepatovac_area ~ s(lipids_PE_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_lipids_PE)
# lipids_PE_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_lipids_PE <- lm(log_prop_hepatovac_area ~ poly(lipids_PE_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_lipids_PE)
# lipids_PE_µg_mg does not significantly affect the response in this model.

poly3_lipids_PE <- lm(log_prop_hepatovac_area ~ poly(lipids_PE_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_lipids_PE)
# lipids_PE_µg_mg does not significantly affect the response in this model.

### total_polar_lipids_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_polar_lipids_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("total_polar_lipids_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_polar_lipids_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_total_polar_lipids_µg_mg = mean(total_polar_lipids_µg_mg),
            sd_total_polar_lipids_µg_mg = sd(total_polar_lipids_µg_mg))

aov01 <- aov(total_polar_lipids_µg_mg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_polar_lipids_µg_mg levels between fish sampled at different season
# Fish sampled in summer present a significantly higher level of total polar lipids (µg.mg1) than fish sampled in winter. 
# No significant difference in total polar lipids (µg.mg1) levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_polar_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_polar_lipids_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_total_polar <- lm(log_prop_hepatovac_area ~ total_polar_lipids_µg_mg + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_polar)
# total_polar_lipids_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_polar_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_total_polar <- gam(log_prop_hepatovac_area ~ s(total_polar_lipids_µg_mg) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_polar)
# total_polar_lipids_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_total_polar <- lm(log_prop_hepatovac_area ~ poly(total_polar_lipids_µg_mg, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_total_polar)
# total_polar_lipids_µg_mg does not significantly affect the response in this model.

poly3_total_polar <- lm(log_prop_hepatovac_area ~ poly(total_polar_lipids_µg_mg, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_total_polar)
# total_polar_lipids_µg_mg does not significantly affect the response in this model.

## plot_polar_lipids -------------------------------------------------
ggarrange(plot_lipids_LPC, plot_lipids_CL, ncol = 2, nrow = 1, labels = c("A","B"))

## total_membrane_lipids_µg_mg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_membrane_lipids_µg_mg) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("total_membrane_lipids_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_membrane_lipids_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_total_membrane_lipids_µg_mg = mean(total_membrane_lipids_µg_mg),
            sd_total_membrane_lipids_µg_mg = sd(total_membrane_lipids_µg_mg))

aov01 <- aov(total_membrane_lipids_µg_mg ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_membrane_lipids_µg_mg levels between fishes sampled at different season and between sex.
# Fish sampled in Seine present a significantly lower level of total_membrane_lipids_µg_mg comparatively to fish sampled in Canche.
# Fish sampled in winter present a significantly higher level of total_membrane_lipids_µg_mg than fish sampled in summer.
# Male present a significantly higher level of total_membrane_lipids_µg_mg than female.
# No significant difference in total_membrane_lipids_µg_mg levels between station.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_membrane_lipids_µg_mg, log_prop_hepatovac_area, col = station) +
  geom_point() +
  geom_smooth(method = lm,  se = F) #+ facet_wrap(~ station)

cor.test(x = data_TOXEM$total_membrane_lipids_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# statistically significant negative correlation

lm_total_membrane <- lm(log_prop_hepatovac_area ~ total_membrane_lipids_µg_mg + station * season + poly(age_month, 2), data = data_clean)
summary(lm_total_membrane)
# total_membrane_lipids_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_membrane_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_total_membrane <- gam(log_prop_hepatovac_area ~ s(total_membrane_lipids_µg_mg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_membrane)
# total_membrane_lipids_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_total_membrane <- lm(log_prop_hepatovac_area ~ poly(total_membrane_lipids_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_total_membrane)
# total_membrane_lipids_µg_mg does not significantly affect the response in this model.

poly3_total_membrane <- lm(log_prop_hepatovac_area ~ poly(total_membrane_lipids_µg_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_total_membrane)
# total_membrane_lipids_µg_mg does not significantly affect the response in this model.

# Compare models
vif(gam_total_membrane)
vif(lm_total_membrane)
# vif < 10 -> no collinearity between variables

AIC(lm_total_membrane, gam_total_membrane)
BIC(lm_total_membrane, gam_total_membrane)

# Plot fitted values
plot(data_TOXEM$total_membrane_lipids_µg_mg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$total_membrane_lipids_µg_mg, fitted(gam_total_membrane), col = "red", pch = 2)
points(data_clean$total_membrane_lipids_µg_mg, fitted(lm_total_membrane), col = "blue", pch = 3)

# Check assumptions
plot(lm_total_membrane)
res <- resid(lm_total_membrane)
plot(fitted(lm_total_membrane), res) # homogeneity of the variances
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res)) # normality of the variances

# lm_total_membrane

## total_lipids_µg_mg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_lipids_µg_mg))
#data_clean <- na.omit(data_TOXEM[, c("total_lipids_µg_mg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_lipids_µg_mg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_total_lipids_µg_mg = mean(total_lipids_µg_mg),
            sd_total_lipids_µg_mg = sd(total_lipids_µg_mg))

aov01 <- aov(total_lipids_µg_mg ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_lipids_µg_mg levels between fishes sampled in the different station and at different season.
# Fish sampled in Seine present a significantly higher level of total_lipids_µg_mg than fish sampled in Canche. 
# Fish sampled in winter present a significantly higher level of total_lipids_µg_mg than fish sampled in summer. 
# No significant difference in total_lipids_µg_mg levels between females and males.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_lipids_µg_mg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_total_lipids <- lm(log_prop_hepatovac_area ~ total_lipids_µg_mg + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_lipids)
# total_lipids_µg_mg does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_lipids_µg_mg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_total_lipids <- gam(log_prop_hepatovac_area ~ s(total_lipids_µg_mg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_lipids)
# total_lipids_µg_mg does not significantly affect the response in this model.

# Polynomial
poly2_total_lipids <- lm(log_prop_hepatovac_area ~ poly(total_lipids_µg_mg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_total_lipids)
# total_lipids_µg_mg does not significantly affect the response in this model.

poly3_total_lipids <- lm(log_prop_hepatovac_area ~ poly(total_lipids_µg_mg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_total_lipids)
# total_lipids_µg_mg does not significantly affect the response in this model.

## HAP_data_TOXEM ------------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### naphtalene_liver_ng_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(naphtalene_liver_ng_g))
#data_clean <- na.omit(data_TOXEM[, c("naphtalene_liver_ng_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = naphtalene_liver_ng_g, fill = sex) +
  #geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_naphtalene_liver_ng_g = mean(naphtalene_liver_ng_g))

aov01 <- aov(naphtalene_liver_ng_g ~ station * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in the level of naphtalene_liver_ng_g between fishes from the different station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(naphtalene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$naphtalene_liver_ng_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_naphtalene <- lm(log_prop_hepatovac_area ~ naphtalene_liver_ng_g + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_naphtalene)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(naphtalene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_naphtalene <- gam(log_prop_hepatovac_area ~ s(naphtalene_liver_ng_g), data = data_TOXEM)
summary(gam_naphtalene)

# Polynomial
poly2_naphtalene <- lm(log_prop_hepatovac_area ~ poly(naphtalene_liver_ng_g, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_naphtalene)

poly3_naphtalene <- lm(log_prop_hepatovac_area ~ poly(naphtalene_liver_ng_g, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_naphtalene)

glm_naphtalene <- glm(log_prop_hepatovac_area ~ naphtalene_liver_ng_g, family=gaussian, data=data_TOXEM)
summary(glm_naphtalene)

### dibenzothiophene_liver_ng_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(dibenzothiophene_liver_ng_g))
#data_clean <- na.omit(data_TOXEM[, c("dibenzothiophene_liver_ng_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = dibenzothiophene_liver_ng_g, color = season, shape = sex) +
  geom_point() +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_naphtalene_liver_ng_g = mean(dibenzothiophene_liver_ng_g))

aov01 <- aov(dibenzothiophene_liver_ng_g ~ station * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in the level of dibenzothiophene_liver_ng_g between fishes from the different station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(dibenzothiophene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$dibenzothiophene_liver_ng_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_dibenzothiophene <- lm(log_prop_hepatovac_area ~ dibenzothiophene_liver_ng_g + station + age_month, data = data_TOXEM)
summary(lm_dibenzothiophene)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(dibenzothiophene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_dibenzothiophene <- gam(log_prop_hepatovac_area ~ s(dibenzothiophene_liver_ng_g), data = data_TOXEM)
summary(gam_dibenzothiophene)

# Polynomial
poly2_dibenzothiophene <- lm(log_prop_hepatovac_area ~ poly(dibenzothiophene_liver_ng_g, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_dibenzothiophene)

poly3_dibenzothiophene <- lm(log_prop_hepatovac_area ~ poly(dibenzothiophene_liver_ng_g, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_dibenzothiophene)

glm_dibenzothiophene <- glm(log_prop_hepatovac_area ~ dibenzothiophene_liver_ng_g, family=gaussian, data=data_clean)

### fluoranthene_liver_ng_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(fluoranthene_liver_ng_g))
#data_clean <- na.omit(data_TOXEM[, c("fluoranthene_liver_ng_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = fluoranthene_liver_ng_g, color = season, shape = sex) +
  geom_point() +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_naphtalene_liver_ng_g = mean(fluoranthene_liver_ng_g))

aov01 <- aov(fluoranthene_liver_ng_g ~ station * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in the level of fluoranthene_liver_ng_g between fishes from the different station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(fluoranthene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$fluoranthene_liver_ng_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_fluoranthene <- lm(log_prop_hepatovac_area ~ fluoranthene_liver_ng_g + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_fluoranthene)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(fluoranthene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_fluoranthene <- gam(log_prop_hepatovac_area ~ s(fluoranthene_liver_ng_g) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_fluoranthene)

# Polynomial
poly2_fluoranthene <- lm(log_prop_hepatovac_area ~ poly(fluoranthene_liver_ng_g, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_fluoranthene)

poly3_fluoranthene <- lm(log_prop_hepatovac_area ~ poly(fluoranthene_liver_ng_g, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_fluoranthene)

glm_fluoranthene <- glm(log_prop_hepatovac_area ~ fluoranthene_liver_ng_g + poly(age_month, 2), family=gaussian, data= data_TOXEM)
summary(glm_fluoranthene)

### pyrene_liver_ng_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(pyrene_liver_ng_g))
#data_clean <- na.omit(data_TOXEM[, c("pyrene_liver_ng_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = pyrene_liver_ng_g, color = season, shape = sex) +
  geom_point() +
  geom_jitter() +
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_pyrene_liver_ng_g = mean(pyrene_liver_ng_g))

aov01 <- aov(pyrene_liver_ng_g ~ station * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in the level of pyrene_liver_ng_g between fishes from the different station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(pyrene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$pyrene_liver_ng_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_pyrene <- lm(log_prop_hepatovac_area ~ pyrene_liver_ng_g + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_pyrene)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(pyrene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_pyrene <- gam(log_prop_hepatovac_area ~ s(pyrene_liver_ng_g) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_pyrene)

# Polynomial
poly2_pyrene <- lm(log_prop_hepatovac_area ~ poly(pyrene_liver_ng_g, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_pyrene)

poly3_pyrene <- lm(log_prop_hepatovac_area ~ poly(pyrene_liver_ng_g, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_pyrene)

glm_pyrene <- glm(log_prop_hepatovac_area ~ pyrene_liver_ng_g + station * season * sex + poly(age_month, 2), family=gaussian, data=data_clean)
summary(glm_pyrene)

### chrysene_liver_ng_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(chrysene_liver_ng_g))
#data_clean <- na.omit(data_TOXEM[, c("chrysene_liver_ng_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = chrysene_liver_ng_g, color = season, shape = sex) +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_chrysene_liver_ng_g = mean(chrysene_liver_ng_g))

aov01 <- aov(chrysene_liver_ng_g ~ station * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in the level of chrysene_liver_ng_g between fishes from the different station.
# Fish sampled in Canche present a significantly higher level of chrysene_liver_ng_g than fish sample in Seine.
# No significant difference in the level of chrysene_liver_ng_g between fishes from the different season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(chrysene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$chrysene_liver_ng_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_chrysene <- lm(log_prop_hepatovac_area ~ chrysene_liver_ng_g + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_chrysene)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(chrysene_liver_ng_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_chrysene <- gam(log_prop_hepatovac_area ~ s(chrysene_liver_ng_g) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_chrysene)

# Polynomial
poly2_chrysene <- lm(log_prop_hepatovac_area ~ poly(chrysene_liver_ng_g, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_chrysene)

poly3_chrysene <- lm(log_prop_hepatovac_area ~ poly(chrysene_liver_ng_g, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_chrysene)

glm_chrysene <- glm(log_prop_hepatovac_area ~ chrysene_liver_ng_g + station + poly(age_month, 2), family=gaussian, data=data_clean)
summary(glm_chrysene)

### total_HAP_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_HAP_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("total_HAP_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_HAP_liver_mg_kg, fill = sex) +
  geom_boxplot() +
  #geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(sex) %>%
  summarise(mean_total_HAP_liver_mg_kg = mean(total_HAP_liver_mg_kg),
            sd_total_HAP_liver_mg_kg = sd(total_HAP_liver_mg_kg))

aov01 <- aov(total_HAP_liver_mg_kg ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_HAP_liver_mg_kg levels between females and males.
# Males present a significantly higher level of total_HAP_liver_mg_kg comparatively to females.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_HAP_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_HAP_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_total_HAP_liver_mg_kg <- lm(log_prop_hepatovac_area ~ total_HAP_liver_mg_kg + sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_HAP_liver_mg_kg)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_HAP_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_total_HAP_liver_mg_kg <- gam(log_prop_hepatovac_area ~ s(total_HAP_liver_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_HAP_liver_mg_kg)

# Polynomial
poly2_total_HAP_liver_mg_kg <- lm(log_prop_hepatovac_area ~ poly(total_HAP_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_total_HAP_liver_mg_kg)

poly3_total_HAP_liver_mg_kg <- lm(log_prop_hepatovac_area ~ poly(total_HAP_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_total_HAP_liver_mg_kg)

glm_total_HAP_liver_mg_kg <- glm(log_prop_hepatovac_area ~ total_HAP_liver_mg_kg + station + poly(age_month, 2), family=gaussian, data=data_clean)
summary(glm_total_HAP_liver_mg_kg)

## PCB_data_TOXEM ------------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### PCB_18 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_18))
#data_clean <- na.omit(data_TOXEM[, c("PCB_18", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_18, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_18 = mean(PCB_18),
            sd_PCB_18 = sd(PCB_18))

aov01 <- aov(PCB_18 ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_18 levels between fishes sampled in the different station.
# Fish sampled in Canche present a significantly higher level of PCB_18 than fish sampled in Seine
# No significant difference in PCB_18 levels between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_18, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_18, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_18 <- lm(log_prop_hepatovac_area ~ PCB_18 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_18)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_18, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_18 <- gam(log_prop_hepatovac_area ~ s(PCB_18), data = data_TOXEM)
summary(gam_PCB_18)

# Polynomial
poly2_PCB_18 <- lm(log_prop_hepatovac_area ~ poly(PCB_18, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_18)

poly3_PCB_18 <- lm(log_prop_hepatovac_area ~ poly(PCB_18, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_18)

### PCB_28 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_28))
#data_clean <- na.omit(data_TOXEM[, c("PCB_28", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_28, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_28 = mean(PCB_28))

aov01 <- aov(PCB_28 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_28 levels between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_28, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_28, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_28 <- lm(log_prop_hepatovac_area ~ PCB_28 + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_28)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_28, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_28 <- gam(log_prop_hepatovac_area ~ s(PCB_28) + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_28)

# Polynomial
poly2_PCB_28 <- lm(log_prop_hepatovac_area ~ poly(PCB_28, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_28)

poly3_PCB_28 <- lm(log_prop_hepatovac_area ~ poly(PCB_28, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_28)

### PCB_31 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_31))
#data_clean <- na.omit(data_TOXEM[, c("PCB_31", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_31, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_PCB_31 = mean(PCB_31),
            sd_PCB_31 = sd(PCB_31))

aov01 <- aov(PCB_31 ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_31 levels between fishes sampled at the different season.
# Fish sampled in winter present a significantly higher level of PCB_31 than fish sampled in summer.
# No significant difference in PCB_31 levels between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_31, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_31, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_31 <- lm(log_prop_hepatovac_area ~ PCB_31 + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_31)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_31, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_31 <- gam(log_prop_hepatovac_area ~ s(PCB_31) + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_31)

# Polynomial
poly2_PCB_31 <- lm(log_prop_hepatovac_area ~ poly(PCB_31, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_31)

poly3_PCB_31 <- lm(log_prop_hepatovac_area ~ poly(PCB_31, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_31)

### PCB_44 *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_44) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("PCB_44", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_44, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_44 = mean(PCB_44))

aov01 <- aov(PCB_44 ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_44 levels between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_44, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_44, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_44 <- lm(log_prop_hepatovac_area ~ PCB_44 + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_44)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_44, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_44 <- gam(log_prop_hepatovac_area ~ s(PCB_44) + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_44)

# Polynomial
poly2_PCB_44 <- lm(log_prop_hepatovac_area ~ poly(PCB_44, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_44)

poly3_PCB_44 <- lm(log_prop_hepatovac_area ~ poly(PCB_44, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_44)

(plot_PCB_44 <- ggplot(data_TOXEM) +
  aes(PCB_44, log_prop_hepatovac_area) + #, col = station, shape = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F))

# Compare models
vif(poly3_PCB_44)

summary(poly3_PCB_44)$r.squared

AIC(poly3_PCB_44)

# Plot fitted values
plot(data_TOXEM$PCB_44, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$PCB_44, fitted(poly3_PCB_44), col = "red", pch = 2)

# Check assumptions
plot(poly3_PCB_44)
res <- resid(poly3_PCB_44)
plot(fitted(poly3_PCB_44), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

# poly3_PCB_44

### PCB_49 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_49))
#data_clean <- na.omit(data_TOXEM[, c("PCB_49", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_49, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_49 = mean(PCB_49),
            sd_PCB_49 = sd(PCB_49))

aov01 <- aov(PCB_49 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_49 levels between station.
# Fish sampled in Seine present a significantly higher level of PCB_49 than fish sampled in Canche.
# No significant difference in PCB_49 levels between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_49, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_49, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_49 <- lm(log_prop_hepatovac_area ~ PCB_49 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_49)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_49, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_49 <- gam(log_prop_hepatovac_area ~ s(PCB_49), data = data_TOXEM)
summary(gam_PCB_49)

# Polynomial
poly2_PCB_49 <- lm(log_prop_hepatovac_area ~ poly(PCB_49, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_49)

poly3_PCB_49 <- lm(log_prop_hepatovac_area ~ poly(PCB_49, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_49)

### PCB_52 *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_52))
#data_clean <- na.omit(data_TOXEM[, c("PCB_52", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_52, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_52 = mean(PCB_52),
            sd_PCB_52 = sd(PCB_52))

aov01 <- aov(PCB_52 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_52 levels between station.
# Fish sampled in Seine present a significantly higher level of PCB_52 than fish sampled in Canche.
# No significant difference in PCB_52 levels between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_52, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_52, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_52 <- lm(log_prop_hepatovac_area ~ PCB_52 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_52)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_52, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_52 <- gam(log_prop_hepatovac_area ~ s(PCB_52) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_52)

# Polynomial
poly2_PCB_52 <- lm(log_prop_hepatovac_area ~ poly(PCB_52, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_52)

poly3_PCB_52 <- lm(log_prop_hepatovac_area ~ poly(PCB_52, 3) + station, data = data_clean)
summary(poly3_PCB_52)

(plot_PCB_52 <-ggplot(data_TOXEM) +
  aes(PCB_52, log_prop_hepatovac_area) + #, col = station, shape = season, group = station) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F))

# Compare models
vif(gam_PCB_52)
vif(poly3_PCB_52)

summary(gam_PCB_52)$r.sq
summary(poly3_PCB_52)$r.squared

AIC(gam_PCB_52, poly3_PCB_52)

# Plot fitted values
plot(data_TOXEM$PCB_52, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$PCB_52, fitted(gam_PCB_52), col = "red", pch = 2)

# Check assumptions
plot(gam_PCB_52)
res <- resid(gam_PCB_52)
plot(fitted(gam_PCB_52), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

# gam_PCB_52

### PCB_77 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_77))
#data_clean <- na.omit(data_TOXEM[, c("PCB_77", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_77, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_77 = mean(PCB_77))

aov01 <- aov(PCB_77 ~ season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")

t.test(PCB_77 ~ season, data = data_TOXEM)
# No significant difference in PCB_77 levels between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_77, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_77, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_77 <- lm(log_prop_hepatovac_area ~ PCB_77 + season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_77)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_77, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_77 <- gam(log_prop_hepatovac_area ~ s(PCB_77), data = data_TOXEM)
summary(gam_PCB_77)

# Polynomial
poly2_PCB_77 <- lm(log_prop_hepatovac_area ~ poly(PCB_77, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_77)

poly3_PCB_77 <- lm(log_prop_hepatovac_area ~ poly(PCB_77, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_77)

### PCB_101 --------------------------------------------------------------
data_clean <- na.omit(data_TOXEM[, c("PCB_101", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_101, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_101 = mean(PCB_101),
            sd_PCB_101 = sd(PCB_101))

aov01 <- aov(PCB_101 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_101 levels between station.
# Fish sampled in Seine present a significantly higher level of PCB_101 than fish sampled in Canche.
# No significant difference in PCB_101 levels between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_101, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_101, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_101 <- lm(log_prop_hepatovac_area ~ PCB_101 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_101)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_101, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_101 <- gam(log_prop_hepatovac_area ~ s(PCB_101) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_101)

# Polynomial
poly2_PCB_101 <- lm(log_prop_hepatovac_area ~ poly(PCB_101, 2) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_101)

poly3_PCB_101 <- lm(log_prop_hepatovac_area ~ poly(PCB_101, 3) + station * season * sex + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_101)

### PCB_105 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_105))
#data_clean <- na.omit(data_TOXEM[, c("PCB_105", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_105, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_PCB_105 = mean(PCB_105),
            sd_PCB_105 = sd(PCB_105))

aov01 <- aov(PCB_105 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_105 levels between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_105, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_105, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_105 <- lm(log_prop_hepatovac_area ~ PCB_105 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_105)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_105, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_105 <- gam(log_prop_hepatovac_area ~ s(PCB_105) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_105)

# Polynomial
poly2_PCB_105 <- lm(log_prop_hepatovac_area ~ poly(PCB_105, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_105)

poly3_PCB_105 <- lm(log_prop_hepatovac_area ~ poly(PCB_105, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_105)

### PCB_110 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_110))
#data_clean <- na.omit(data_TOXEM[, c("PCB_110", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_110, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_110 = mean(PCB_110),
            sd_PCB_110 = sd(PCB_110))

aov01 <- aov(PCB_110 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_110 levels between station.
# Fish sampled in Seine present a significantly higher PCB_110 level compared to fish sampled in Canche.
# No significant difference in PCB_110 level between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_110, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_110, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_110 <- lm(log_prop_hepatovac_area ~ PCB_110 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_110)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_110, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_110 <- gam(log_prop_hepatovac_area ~ s(PCB_110) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_110)

# Polynomial
poly2_PCB_110 <- lm(log_prop_hepatovac_area ~ poly(PCB_110, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_110)

poly3_PCB_110 <- lm(log_prop_hepatovac_area ~ poly(PCB_110, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_110)

### PCB_118 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_118))
#data_clean <- na.omit(data_TOXEM[, c("PCB_118", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_118, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_118 = mean(PCB_118),
            sd_PCB_118 = sd(PCB_118))

aov01 <- aov(PCB_118 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_118 levels between station.
# Fish sampled in Seine present a significantly higher PCB_118 level compared to fish sampled in Canche.
# No significant difference in PCB_118 level between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_118, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_118, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_118 <- lm(log_prop_hepatovac_area ~ PCB_118 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_118)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_118, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_118 <- gam(log_prop_hepatovac_area ~ s(PCB_118) + station * season * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_118)

# Polynomial
poly2_PCB_118 <- lm(log_prop_hepatovac_area ~ poly(PCB_118, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_118)

poly3_PCB_118 <- lm(log_prop_hepatovac_area ~ poly(PCB_118, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_118)

### PCB_128 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_128))
#data_clean <- na.omit(data_TOXEM[, c("PCB_128", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_128, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_128 = mean(PCB_128))

aov01 <- aov(PCB_128 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_128 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_128, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_128, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_128 <- lm(log_prop_hepatovac_area ~ PCB_128 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_128)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_128, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_128 <- gam(log_prop_hepatovac_area ~ s(PCB_128) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_128)

# Polynomial
poly2_PCB_128 <- lm(log_prop_hepatovac_area ~ poly(PCB_128, 2) + station * season, data = data_clean)
summary(poly2_PCB_128)

poly3_PCB_128 <- lm(log_prop_hepatovac_area ~ poly(PCB_128, 3) + station * season, data = data_clean)
summary(poly3_PCB_128)

### PCB_132 *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_132))
#data_clean <- na.omit(data_TOXEM[, c("PCB_132", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_132, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_132 = mean(PCB_132))

aov01 <- aov(PCB_132 ~ season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_132 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_132, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_132, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_132 <- lm(log_prop_hepatovac_area ~ PCB_132 + season * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_132)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_132, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_132 <- gam(log_prop_hepatovac_area ~ s(PCB_132) + season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_132)

# Polynomial
poly2_PCB_132 <- lm(log_prop_hepatovac_area ~ poly(PCB_132, 2) + season * sex + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_132)

(plot_PCB_132 <-ggplot(data_TOXEM) +
  aes(PCB_132, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F))

poly3_PCB_132 <- lm(log_prop_hepatovac_area ~ poly(PCB_132, 3) + season * sex+ poly(age_month, 2), data = data_clean)
summary(poly3_PCB_132)

# Compare models
vif(poly2_PCB_132)

summary(poly2_PCB_132)$r.squared

AIC(poly2_PCB_132)

# Plot fitted values
plot(data_TOXEM$PCB_132, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$PCB_132, fitted(poly2_PCB_132), col = "red", pch = 2)

# Check assumptions
plot(poly2_PCB_132)
res <- resid(poly2_PCB_132)
plot(fitted(poly2_PCB_132), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

### PCB_135 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_135))
#data_clean <- na.omit(data_TOXEM[, c("PCB_135", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_135, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_135 = mean(PCB_135))

aov01 <- aov(PCB_135 ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_135 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_135, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_135, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_135 <- lm(log_prop_hepatovac_area ~ PCB_135 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_135)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_135, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_135 <- gam(log_prop_hepatovac_area ~ s(PCB_135) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_135)

# Polynomial
poly2_PCB_135 <- lm(log_prop_hepatovac_area ~ poly(PCB_135, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_135)

poly3_PCB_135 <- lm(log_prop_hepatovac_area ~ poly(PCB_135, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_135)

### PCB_138 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_138))
#data_clean <- na.omit(data_TOXEM[, c("PCB_138", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_138, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_138 = mean(PCB_138))

aov01 <- aov(PCB_138 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_138 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_138, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_138, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_138 <- lm(log_prop_hepatovac_area ~ PCB_138 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_138)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_138, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_138 <- gam(log_prop_hepatovac_area ~ s(PCB_138) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_138)

# Polynomial
poly2_PCB_138 <- lm(log_prop_hepatovac_area ~ poly(PCB_138, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_138)

poly3_PCB_138 <- lm(log_prop_hepatovac_area ~ poly(PCB_138, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_138)

### PCB_149 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_149))
#data_clean <- na.omit(data_TOXEM[, c("PCB_149", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_149, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_149 = mean(PCB_149),
            sd_PCB_149 = sd(PCB_149))

aov01 <- aov(PCB_149 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_149 level between station.
# Fish sampled in Seine present significantly higher level of PCB_149 comparatively to fish sampled in Canche.
# No significant difference in PCB_149 level between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_149, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_149, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_149 <- lm(log_prop_hepatovac_area ~ PCB_149 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_149)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_149, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_149 <- gam(log_prop_hepatovac_area ~ s(PCB_149) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_149)

# Polynomial
poly2_PCB_149 <- lm(log_prop_hepatovac_area ~ poly(PCB_149, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_149)

poly3_PCB_149 <- lm(log_prop_hepatovac_area ~ poly(PCB_149, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_149)

### PCB_153 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_153))
#data_clean <- na.omit(data_TOXEM[, c("PCB_153", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_153, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_PCB_153 = mean(PCB_153),
            sd_PCB_153 = sd(PCB_153))

aov01 <- aov(PCB_153 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in PCB_153 level between station.
# Fish sampled in Seine present significantly higher level of PCB_153 comparatively to fish sampled in Canche.
# No significant difference in PCB_153 level between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_153, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_153, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_153 <- lm(log_prop_hepatovac_area ~ PCB_153 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_153)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_153, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_153 <- gam(log_prop_hepatovac_area ~ s(PCB_153) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_153)

# Polynomial
poly2_PCB_153 <- lm(log_prop_hepatovac_area ~ poly(PCB_153, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_153)

poly3_PCB_153 <- lm(log_prop_hepatovac_area ~ poly(PCB_153, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_153)

### PCB_156 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_156))
#data_clean <- na.omit(data_TOXEM[, c("PCB_156", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_156, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_156 = mean(PCB_156))

aov01 <- aov(PCB_156 ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_156 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_156, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_156, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_156 <- lm(log_prop_hepatovac_area ~ PCB_156 + season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_156)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_156, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_156 <- gam(log_prop_hepatovac_area ~ s(PCB_156), data = data_TOXEM)
summary(gam_PCB_156)

# Polynomial
poly2_PCB_156 <- lm(log_prop_hepatovac_area ~ poly(PCB_156, 2), data = data_clean)
summary(poly2_PCB_156)

poly3_PCB_156 <- lm(log_prop_hepatovac_area ~ poly(PCB_156, 3), data = data_clean)
summary(poly3_PCB_156)

### PCB_169 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_169))
#data_clean <- na.omit(data_TOXEM[, c("PCB_169", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_169, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_169 = mean(PCB_169))

aov01 <- aov(PCB_169 ~ sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_169 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_169, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_169, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_169 <- lm(log_prop_hepatovac_area ~ PCB_169, data = data_TOXEM)
summary(lm_PCB_169)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_169, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_169 <- gam(log_prop_hepatovac_area ~ s(PCB_169), data = data_TOXEM)
summary(gam_PCB_169)

# Polynomial
poly2_PCB_169 <- lm(log_prop_hepatovac_area ~ poly(PCB_169, 2), data = data_clean)
summary(poly2_PCB_169)

poly3_PCB_169 <- lm(log_prop_hepatovac_area ~ poly(PCB_169, 3), data = data_clean)
summary(poly3_PCB_169)

### PCB_170 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_170))
#data_clean <- na.omit(data_TOXEM[, c("PCB_170", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_170, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_170 = mean(PCB_170))

aov01 <- aov(PCB_170 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_170 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_170, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_170, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_170 <- lm(log_prop_hepatovac_area ~ PCB_170 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_170)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_170, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_170 <- gam(log_prop_hepatovac_area ~ s(PCB_170) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_170)

# Polynomial
poly2_PCB_170 <- lm(log_prop_hepatovac_area ~ poly(PCB_170, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_170)

poly3_PCB_170 <- lm(log_prop_hepatovac_area ~ poly(PCB_170, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_170)

### PCB_180 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_180))
#data_clean <- na.omit(data_TOXEM[, c("PCB_180", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_180, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_180 = mean(PCB_180))

aov01 <- aov(PCB_180 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_180 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_180, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_180, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_180 <- lm(log_prop_hepatovac_area ~ PCB_180 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_180)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_180, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_180 <- gam(log_prop_hepatovac_area ~ s(PCB_180) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_180)

# Polynomial
poly2_PCB_180 <- lm(log_prop_hepatovac_area ~ poly(PCB_180, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_180)

poly3_PCB_180 <- lm(log_prop_hepatovac_area ~ poly(PCB_180, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_180)

### PCB_187 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(PCB_187))
#data_clean <- na.omit(data_TOXEM[, c("PCB_187", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = PCB_187, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_PCB_187 = mean(PCB_187))

aov01 <- aov(PCB_187 ~ station * season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in PCB_187 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(PCB_187, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$PCB_187, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_PCB_187 <- lm(log_prop_hepatovac_area ~ PCB_187 + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_PCB_187)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(PCB_187, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_PCB_187 <- gam(log_prop_hepatovac_area ~ s(PCB_187) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_PCB_187)

# Polynomial
poly2_PCB_187 <- lm(log_prop_hepatovac_area ~ poly(PCB_187, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_PCB_187)

poly3_PCB_187 <- lm(log_prop_hepatovac_area ~ poly(PCB_187, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_PCB_187)

### total_PCB_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_PCB_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("total_PCB_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_PCB_liver_mg_kg, fill = sex) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(sex) %>%
  summarise(mean_total_PCB_liver_mg_kg = mean(total_PCB_liver_mg_kg),
            sd_total_PCB_liver_mg_kg = sd(total_PCB_liver_mg_kg))

aov01 <- aov(total_PCB_liver_mg_kg ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_PCB_liver_mg_kg level measured in fish liver between the two station.
# Fish from Seine present significantly higher total_PCB_liver_mg_kg levels comparatively to fish sampled in Canche.  
# No significant difference in total_PCB_liver_mg_kg level between season and sex.

# Les analyses de PCB et de pesticides font clairement ressortir des concentrations environ 10 fois plus fortes en PCB dans le foie en Seine vs Canche, quelle que soit la saison ; les teneurs en pesticides étant 40 fois et 4 fois plus fortes, respectivement en hiver et en été, en Seine vs Canche.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_PCB_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$total_PCB_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_total_PCB <- lm(log_prop_hepatovac_area ~ total_PCB_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_PCB)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_PCB_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_total_PCB <- gam(log_prop_hepatovac_area ~ s(total_PCB_liver_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_PCB)

# Polynomial
poly2_total_PCB <- lm(log_prop_hepatovac_area ~ poly(total_PCB_liver_mg_kg, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_total_PCB)

poly3_total_PCB <- lm(log_prop_hepatovac_area ~ poly(total_PCB_liver_mg_kg, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_total_PCB)

## plot_PCB -------------------------------------------------
ggarrange(plot_PCB_44, plot_PCB_52, plot_PCB_132, ncol = 2, nrow = 2, labels = c("A","B","C"))

## pesticides_data_TOXEM ------------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### four_four_dde ***--------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(four_four_dde) & !is.na(age_month))
#data_clean <- na.omit(data_TOXEM[, c("four_four_dde", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = four_four_dde, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_four_four_dde = mean(four_four_dde),
            sd_four_four_dde = sd(four_four_dde))

aov01 <- aov(four_four_dde ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in four_four_dde level measured in fish liver between the two station.
# Fish from Seine present significantly higher four_four_dde levels comparatively to fish sampled in Canche.  
# No significant difference in four_four_dde level between season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(four_four_dde, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$four_four_dde, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_four_four_dde <- lm(log_prop_hepatovac_area ~ four_four_dde + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_four_four_dde)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(four_four_dde, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = T)

gam_four_four_dde <- gam(log_prop_hepatovac_area ~ s(four_four_dde) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_four_four_dde)

# Polynomial
poly2_four_four_dde <- lm(log_prop_hepatovac_area ~ poly(four_four_dde, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_four_four_dde)

poly3_four_four_dde <- lm(log_prop_hepatovac_area ~ poly(four_four_dde, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_four_four_dde)

ggplot(data_TOXEM) +
  aes(four_four_dde, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = F)

# Compare models
vif(gam_four_four_dde)
vif(poly2_four_four_dde)
vif(poly3_four_four_dde)

summary(gam_four_four_dde)$r.sq
summary(poly2_four_four_dde)$r.squared
summary(poly3_four_four_dde)$r.squared

AIC(gam_four_four_dde, poly2_four_four_dde, poly3_four_four_dde)

anova(gam_four_four_dde, poly2_four_four_dde, poly3_four_four_dde, test = "F")

# Plot fitted values
plot(data_TOXEM$four_four_dde, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$four_four_dde, fitted(gam_four_four_dde), col = "green", pch = 4)
points(data_clean$four_four_dde, fitted(poly3_four_four_dde), col = "blue", pch = 3)
points(data_clean$four_four_dde, fitted(poly2_four_four_dde), col = "red", pch = 2)

# Check assumptions
plot(poly3_four_four_dde)
res <- resid(poly3_four_four_dde)
plot(fitted(poly3_four_four_dde), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

# poly3_four_four_dde

### four_four_ddd --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(four_four_ddd))
#data_clean <- na.omit(data_TOXEM[, c("four_four_ddd", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = four_four_ddd, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_four_four_ddd = mean(four_four_ddd))

aov01 <- aov(four_four_ddd ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in four_four_ddd level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(four_four_ddd, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$four_four_ddd, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_four_four_ddd <- lm(log_prop_hepatovac_area ~ four_four_ddd + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_four_four_ddd)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(four_four_ddd, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_four_four_ddd <- gam(log_prop_hepatovac_area ~ s(four_four_ddd), data = data_TOXEM)
summary(gam_four_four_ddd)

# Polynomial
poly2_four_four_ddd <- lm(log_prop_hepatovac_area ~ poly(four_four_ddd, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_four_four_ddd)

poly3_four_four_ddd <- lm(log_prop_hepatovac_area ~ poly(four_four_ddd, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_four_four_ddd)

### total_pesticides_liver_mg_kg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_pesticides_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("total_pesticides_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_pesticides_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_total_pesticides_liver_mg_kg = mean(total_pesticides_liver_mg_kg),
            sd_total_pesticides_liver_mg_kg = sd(total_pesticides_liver_mg_kg))

aov01 <- aov(total_pesticides_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_pesticides_liver_mg_kg level between station and season.
# Fish sampled in Seine present a significantly higher level of total_pesticides_liver_mg_kg comparatively to fish sampled in Canche.
# Fish sampled in winter present a significantly higher level of total_pesticides_liver_mg_kg comparatively to fish sampled in summer.
# No significant difference in total_pesticides_liver_mg_kg level between sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_pesticides_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_pesticides_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_total_pesticides <- lm(log_prop_hepatovac_area ~ total_pesticides_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_pesticides)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_pesticides_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~ station)

gam_total_pesticides <- gam(log_prop_hepatovac_area ~ s(total_pesticides_liver_mg_kg) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_total_pesticides)

# Polynomial
poly2_total_pesticides <- lm(log_prop_hepatovac_area ~ poly(total_pesticides_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_total_pesticides)

poly3_total_pesticides <- lm(log_prop_hepatovac_area ~ poly(total_pesticides_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_total_pesticides)

## PBDE (polybromodiphenylether) -------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### BDE_47 --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(BDE_47))
#data_clean <- na.omit(data_TOXEM[, c("BDE_47", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = BDE_47, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_BDE_47 = mean(BDE_47))

aov01 <- aov(BDE_47 ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in BDE_47 level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(BDE_47, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$BDE_47, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_BDE_47 <- lm(log_prop_hepatovac_area ~ BDE_47 + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_BDE_47)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(BDE_47, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F) +
  facet_wrap(~ station)

gam_BDE_47 <- gam(log_prop_hepatovac_area ~ s(BDE_47), data = data_TOXEM)
summary(gam_BDE_47)

# Polynomial
poly2_BDE_47 <- lm(log_prop_hepatovac_area ~ poly(BDE_47, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_BDE_47)

poly3_BDE_47 <- lm(log_prop_hepatovac_area ~ poly(BDE_47, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_BDE_47)

### total_PBDE_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_PBDE_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("total_PBDE_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_PBDE_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(season) %>%
  summarise(mean_total_PBDE_liver_mg_kg = mean(total_PBDE_liver_mg_kg),
            sd_total_PBDE_liver_mg_kg = sd(total_PBDE_liver_mg_kg))

aov01 <- aov(total_PBDE_liver_mg_kg ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in total_PBDE_liver_mg_kg level between season.
# Fish sampled in summer present a significantly higher level of total_PBDE_liver_mg_kg measured in liver comparatively to fish sampled in winter.
# No significant difference in total_PBDE_liver_mg_kg level between station and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_PBDE_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) + 
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_PBDE_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_total_PBDE <- lm(log_prop_hepatovac_area ~ total_PBDE_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_total_PBDE)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_PBDE_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_total_PBDE <- gam(log_prop_hepatovac_area ~ s(total_PBDE_liver_mg_kg), data = data_TOXEM)
summary(gam_total_PBDE)

# Polynomial
poly2_total_PBDE <- lm(log_prop_hepatovac_area ~ poly(total_PBDE_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_total_PBDE)

poly3_total_PBDE <- lm(log_prop_hepatovac_area ~ poly(total_PBDE_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_total_PBDE)

### hydroxypyren_µg_g --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(hydroxypyren_µg_g))
#data_clean <- na.omit(data_TOXEM[, c("hydroxypyren_µg_g", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = hydroxypyren_µg_g, fill = sex) +
  geom_boxplot() +
  #geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season) %>%
  summarise(mean_hydroxypyren_µg_g = mean(hydroxypyren_µg_g),
            sd_hydroxypyren_µg_g = sd(hydroxypyren_µg_g))

aov01 <- aov(hydroxypyren_µg_g ~ station * season * sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in the level of hydroxypyren_µg_g between fishes from the different station and season.
# Fish sampled in Seine present a significantly higher level of hydroxypyren_µg_g than fish sample in Canche.
# Fish sampled in winter present a significantly higher level of hydroxypyren_µg_g than fish sample in summer.
# No significant difference in the level of hydroxypyren_µg_g between females and males.

# Linear regression
ggplot(data_TOXEM) +
  aes(hydroxypyren_µg_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   facet_wrap(~ station)

cor.test(x = data_TOXEM$hydroxypyren_µg_g, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")
# no statistically significant correlation

lm_hydroxypyren <- lm(log_prop_hepatovac_area ~ hydroxypyren_µg_g + station * season + poly(age_month, 2), data = data_TOXEM)
summary(lm_hydroxypyren)
# hydroxypyren_µg_g does not significantly affect the response in this model.

# Non linear regression
# GAM
ggplot(data_TOXEM) +
  aes(hydroxypyren_µg_g, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = NULL, se = F) + facet_wrap(~ station)

gam_hydroxypyren <- gam(log_prop_hepatovac_area ~ s(hydroxypyren_µg_g) + station * season + poly(age_month, 2), data = data_TOXEM)
summary(gam_hydroxypyren)
# hydroxypyren_µg_g does not significantly affect the response in this model.

# Polynomial
poly2_hydroxypyren <- lm(log_prop_hepatovac_area ~ poly(hydroxypyren_µg_g, 2) + station * season + poly(age_month, 2), data = data_clean)
summary(poly2_hydroxypyren)
# hydroxypyren_µg_g does not significantly affect the response in this model.

poly3_hydroxypyren <- lm(log_prop_hepatovac_area ~ poly(hydroxypyren_µg_g, 3) + station * season + poly(age_month, 2), data = data_clean)
summary(poly3_hydroxypyren)
# hydroxypyren_µg_g does not significantly affect the response in this model.

glm_hydroxypyren <- glm(log_prop_hepatovac_area ~ hydroxypyren_µg_g + station * season * poly(age_month, 2), family=gaussian, data=data_clean)
summary(glm_hydroxypyren)

## metals_data_TOXEM ------------------------------------------------------
rm(list = ls()[! ls() %in% c("data_TOXEM")])

### Ag_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Ag_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Ag_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Ag_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Ag_liver_mg_kg = mean(Ag_liver_mg_kg))

aov01 <- aov(Ag_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Ag_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Ag_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Ag_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Ag <- lm(log_prop_hepatovac_area ~ Ag_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Ag)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Ag_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Ag <- gam(log_prop_hepatovac_area ~ s(Ag_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Ag)

# Polynomial
poly2_Ag <- lm(log_prop_hepatovac_area ~ poly(Ag_liver_mg_kg, 2) + station * sex, data = data_clean)
summary(poly2_Ag)

poly3_Ag <- lm(log_prop_hepatovac_area ~ poly(Ag_liver_mg_kg, 3) + station * sex, data = data_clean)
summary(poly3_Ag)

### Al_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Al_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Al_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Al_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Al_liver_mg_kg = mean(Al_liver_mg_kg))

aov01 <- aov(Al_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Al_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Al_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Al_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Al <- lm(log_prop_hepatovac_area ~ Al_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Al)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Al_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Al <- gam(log_prop_hepatovac_area ~ s(Al_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Al)

# Polynomial
poly2_Al <- lm(log_prop_hepatovac_area ~ poly(Al_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Al)

poly3_Al <- lm(log_prop_hepatovac_area ~ poly(Al_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Al)

### Cd_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Cd_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Cd_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Cd_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Cd_liver_mg_kg = mean(Cd_liver_mg_kg))

aov01 <- aov(Cd_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Cd_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Cd_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Cd_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Cd <- lm(log_prop_hepatovac_area ~ Cd_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Cd)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Cd_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Cd <- gam(log_prop_hepatovac_area ~ s(Cd_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Cd)

# Polynomial
poly2_Cd <- lm(log_prop_hepatovac_area ~ poly(Cd_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Cd)

poly3_Cd <- lm(log_prop_hepatovac_area ~ poly(Cd_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Cd)

### Co_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Co_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Co_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Co_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Co_liver_mg_kg = mean(Co_liver_mg_kg))

aov01 <- aov(Co_liver_mg_kg ~ sex, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Co_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Co_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Co_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Co <- lm(log_prop_hepatovac_area ~ Co_liver_mg_kg + sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Co)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Co_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Co <- gam(log_prop_hepatovac_area ~ s(Co_liver_mg_kg), data = data_TOXEM)
summary(gam_Co)

# Polynomial
poly2_Co <- lm(log_prop_hepatovac_area ~ poly(Co_liver_mg_kg, 2) + sex + poly(age_month, 2), data = data_clean)
summary(poly2_Co)

poly3_Co <- lm(log_prop_hepatovac_area ~ poly(Co_liver_mg_kg, 3) + sex + poly(age_month, 2), data = data_clean)
summary(poly3_Co)

### Cr_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Cr_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Cr_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Cr_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Cr_liver_mg_kg = mean(Cr_liver_mg_kg))

aov01 <- aov(Cr_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Cr_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Cr_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Cr_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Cr <- lm(log_prop_hepatovac_area ~ Cr_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Cr)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Cr_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Cr <- gam(log_prop_hepatovac_area ~ s(Cr_liver_mg_kg), data = data_TOXEM)
summary(gam_Cr)

# Polynomial
poly2_Cr <- lm(log_prop_hepatovac_area ~ poly(Cr_liver_mg_kg, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_Cr)

poly3_Cr <- lm(log_prop_hepatovac_area ~ poly(Cr_liver_mg_kg, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_Cr)

### Cu_liver_mg_kg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Cu_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Cu_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Cu_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Cu_liver_mg_kg = mean(Cu_liver_mg_kg))

aov01 <- aov(Cu_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Cu_liver_mg_kg level between station, season and sex.

# Linear regression
(Cu_liver <- ggplot(data_TOXEM) +
  aes(Cu_liver_mg_kg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = F))

cor.test(x = data_TOXEM$Cu_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Cu <- lm(log_prop_hepatovac_area ~ Cu_liver_mg_kg, data = data_TOXEM)
summary(lm_Cu)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Cu_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Cu <- gam(log_prop_hepatovac_area ~ s(Cu_liver_mg_kg) + station * sex, data = data_TOXEM)
summary(gam_Cu)

# Polynomial
poly2_Cu <- lm(log_prop_hepatovac_area ~ poly(Cu_liver_mg_kg, 2) + station * sex, data = data_clean)
summary(poly2_Cu)

poly3_Cu <- lm(log_prop_hepatovac_area ~ poly(Cu_liver_mg_kg, 3) + station * sex, data = data_clean)
summary(poly3_Cu)

# Compare models
vif(lm_Cu)

summary(lm_Cu)$r.squared

AIC(lm_Cu)

# Plot fitted values
plot(data_TOXEM$Cu_liver_mg_kg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$Cu_liver_mg_kg, fitted(lm_Cu), col = "green", pch = 4)

# Check assumptions
plot(lm_Cu)
res <- resid(lm_Cu)
plot(fitted(lm_Cu), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

### Fe_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Fe_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Fe_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Fe_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Fe_liver_mg_kg = mean(Fe_liver_mg_kg))

aov01 <- aov(Fe_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Fe_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Fe_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Fe_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Fe <- lm(log_prop_hepatovac_area ~ Fe_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Fe)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Fe_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Fe <- gam(log_prop_hepatovac_area ~ s(Fe_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Fe)

# Polynomial
poly2_Fe <- lm(log_prop_hepatovac_area ~ poly(Fe_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Fe)

poly3_Fe <- lm(log_prop_hepatovac_area ~ poly(Fe_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Fe)

### Hg_liver_mg_kg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Hg_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Hg_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Hg_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station) %>%
  summarise(mean_Hg_liver_mg_kg = mean(Hg_liver_mg_kg),
            sd_Hg_liver_mg_kg = sd(Hg_liver_mg_kg))

aov01 <- aov(Hg_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# Significant difference in Hg_liver_mg_kg level between station.
# Fish from Seine present a significantly higher level of Hg_liver_mg_kg measured in liver comparatively to fish from Canche.
# No significant difference in Hg_liver_mg_kg level between season and sex.

# Linear regression
(Hg_liver <- ggplot(data_TOXEM) +
  aes(Hg_liver_mg_kg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = F))

cor.test(x = data_TOXEM$Hg_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Hg <- lm(log_prop_hepatovac_area ~ Hg_liver_mg_kg + poly(age_month, 2), data = data_TOXEM)
summary(lm_Hg)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Hg_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Hg <- gam(log_prop_hepatovac_area ~ s(Hg_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Hg)

# Polynomial
poly2_Hg <- lm(log_prop_hepatovac_area ~ poly(Hg_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Hg)

poly3_Hg <- lm(log_prop_hepatovac_area ~ poly(Hg_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Hg)

# Compare models
vif(lm_Hg)

summary(lm_Hg)$r.squared

AIC(lm_Hg)

# Plot fitted values
plot(data_TOXEM$Hg_liver_mg_kg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$Hg_liver_mg_kg, fitted(lm_Hg), col = "green", pch = 4)

# Check assumptions
plot(lm_Hg)
res <- resid(lm_Hg)
plot(fitted(lm_Hg), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

### Mg_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Mg_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Mg_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Mg_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Mg_liver_mg_kg = mean(Mg_liver_mg_kg))

aov01 <- aov(Mg_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Mg_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Mg_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Mg_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Mg <- lm(log_prop_hepatovac_area ~ Mg_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Mg)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Mg_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_sMgoth(method = "gam", se = F)

gam_Mg <- gam(log_prop_hepatovac_area ~ s(Mg_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Mg)

# Polynomial
poly2_Mg <- lm(log_prop_hepatovac_area ~ poly(Mg_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Mg)

poly3_Mg <- lm(log_prop_hepatovac_area ~ poly(Mg_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Mg)

### Mn_liver_mg_kg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Mn_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Mn_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Mn_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Mn_liver_mg_kg = mean(Mn_liver_mg_kg))

aov01 <- aov(Mn_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Mn_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Mn_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Mn_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Mn <- lm(log_prop_hepatovac_area ~ Mn_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Mn)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Mn_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Mn <- gam(log_prop_hepatovac_area ~ s(Mn_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Mn)

# Polynomial
poly2_Mn <- lm(log_prop_hepatovac_area ~ poly(Mn_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_Mn)

(Mn_liver <- ggplot(data_TOXEM) +
  aes(Mn_liver_mg_kg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = F))

poly3_Mn <- lm(log_prop_hepatovac_area ~ poly(Mn_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_Mn)

# Compare models
vif(poly2_Mn)
vif(poly3_Mn)

summary(poly2_Mn)$r.squared
summary(poly3_Mn)$r.squared

AIC(poly2_Mn, poly3_Mn)

anova(poly2_Mn, poly3_Mn, test = "F")

# Plot fitted values
plot(data_TOXEM$Mn_liver_mg_kg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_TOXEM$Mn_liver_mg_kg, fitted(poly2_Mn), col = "green", pch = 4)
points(data_TOXEM$Mn_liver_mg_kg, fitted(poly3_Mn), col = "red", pch = 2)

# Check assumptions
plot(poly2_Mn)
res <- resid(poly2_Mn)
plot(fitted(poly2_Mn), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

### Mo_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Mo_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Mo_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Mo_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Mo_liver_mg_kg = mean(Mo_liver_mg_kg))

aov01 <- aov(Mo_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Mo_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Mo_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Mo_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Mo <- lm(log_prop_hepatovac_area ~ Mo_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Mo)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Mo_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Mo <- gam(log_prop_hepatovac_area ~ s(Mo_liver_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_Mo)

# Polynomial
poly2_Mo <- lm(log_prop_hepatovac_area ~ poly(Mo_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_Mo)

poly3_Mo <- lm(log_prop_hepatovac_area ~ poly(Mo_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_Mo)

### Ni_liver_mg_kg *** --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Ni_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Ni_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Ni_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Ni_liver_mg_kg = mean(Ni_liver_mg_kg))

aov01 <- aov(Ni_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Ni_liver_mg_kg level between station, season and sex.

# Linear regression
(Ni_liver <- ggplot(data_TOXEM) +
  aes(Ni_liver_mg_kg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = lm,  se = F))

cor.test(x = data_TOXEM$Ni_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Ni <- lm(log_prop_hepatovac_area ~ Ni_liver_mg_kg + poly(age_month, 2), data = data_TOXEM)
summary(lm_Ni)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Ni_liver_mg_kg, log_prop_hepatovac_area) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Ni <- gam(log_prop_hepatovac_area ~ s(Ni_liver_mg_kg) + poly(age_month, 2), data = data_TOXEM)
summary(gam_Ni)

# Polynomial
poly2_Ni <- lm(log_prop_hepatovac_area ~ poly(Ni_liver_mg_kg, 2) + poly(age_month, 2), data = data_clean)
summary(poly2_Ni)

poly3_Ni <- lm(log_prop_hepatovac_area ~ poly(Ni_liver_mg_kg, 3) + poly(age_month, 2), data = data_clean)
summary(poly3_Ni)

# Compare models
vif(lm_Ni)
vif(gam_Ni)
vif(poly3_Ni)

summary(lm_Ni)$r.squared
summary(gam_Ni)$r.sq
summary(poly3_Ni)$r.squared

AIC(lm_Ni, gam_Ni, poly3_Ni)

anova(lm_Ni, gam_Ni, poly3_Ni, test = "F")

# Plot fitted values
plot(data_TOXEM$Ni_liver_mg_kg, data_TOXEM$log_prop_hepatovac_area, main = "Actual vs Fitted")
points(data_clean$Ni_liver_mg_kg, fitted(lm_Ni), col = "green", pch = 4)
points(data_clean$Ni_liver_mg_kg, fitted(gam_Ni), col = "blue", pch = 3)
points(data_clean$Ni_liver_mg_kg, fitted(poly3_Ni), col = "red", pch = 2)

# Check assumptions
plot(lm_Ni)
res <- resid(lm_Ni)
plot(fitted(lm_Ni), res)
abline(0, 0)
qqnorm(res)
qqline(res)
plot(density(res))

### Pb_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Pb_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Pb_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Pb_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Pb_liver_mg_kg = mean(Pb_liver_mg_kg))

aov01 <- aov(Pb_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Pb_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Pb_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Pb_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Pb <- lm(log_prop_hepatovac_area ~ Pb_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Pb)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Pb_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_sPboth(method = "gam", se = F)

gam_Pb <- gam(log_prop_hepatovac_area ~ s(Pb_liver_mg_kg), data = data_TOXEM)
summary(gam_Pb)

# Polynomial
poly2_Pb <- lm(log_prop_hepatovac_area ~ poly(Pb_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_Pb)

poly3_Pb <- lm(log_prop_hepatovac_area ~ poly(Pb_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_Pb)

### Sr_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Sr_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Sr_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Sr_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Sr_liver_mg_kg = mean(Sr_liver_mg_kg))

aov01 <- aov(Sr_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Sr_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Sr_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Sr_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Sr <- lm(log_prop_hepatovac_area ~ Sr_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Sr)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Sr_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Sr <- gam(log_prop_hepatovac_area ~ s(Sr_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Sr)

# Polynomial
poly2_Sr <- lm(log_prop_hepatovac_area ~ poly(Sr_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Sr)

poly3_Sr <- lm(log_prop_hepatovac_area ~ poly(Sr_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Sr)

### Ti_liver_mg_kg --------------------------------------------------------------
data_clean <- na.omit(data_TOXEM[, c("Ti_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Ti_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Ti_liver_mg_kg = mean(Ti_liver_mg_kg))

aov01 <- aov(Ti_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Ti_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Ti_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Ti_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Ti <- lm(log_prop_hepatovac_area ~ Ti_liver_mg_kg + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(lm_Ti)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Ti_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_sTioth(method = "gam", se = F)

gam_Ti <- gam(log_prop_hepatovac_area ~ s(Ti_liver_mg_kg) + station * sex + poly(age_month, 2), data = data_TOXEM)
summary(gam_Ti)

# Polynomial
ggplot(data_TOXEM) +
  aes(Ti_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_sTioth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "blue") +
  geom_sTioth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "red")

poly2_Ti <- lm(log_prop_hepatovac_area ~ poly(Ti_liver_mg_kg, 2) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly2_Ti)

poly3_Ti <- lm(log_prop_hepatovac_area ~ poly(Ti_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Ti)

### V_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(V_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("V_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = V_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_V_liver_mg_kg = mean(V_liver_mg_kg))

aov01 <- aov(V_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in V_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(V_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +   
  facet_wrap(~ station)

cor.test(x = data_TOXEM$V_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_V <- lm(log_prop_hepatovac_area ~ V_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_V)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(V_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_sVoth(method = "gam", se = F)

gam_V <- gam(log_prop_hepatovac_area ~ s(V_liver_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_V)

# Polynomial
poly2_V <- lm(log_prop_hepatovac_area ~ poly(V_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_V)

poly3_V <- lm(log_prop_hepatovac_area ~ poly(V_liver_mg_kg, 3) + station + poly(age_month, 2), data = data_clean)
summary(poly3_V)

### Zn_liver_mg_kg --------------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(Zn_liver_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("Zn_liver_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = Zn_liver_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_Zn_liver_mg_kg = mean(Zn_liver_mg_kg))

aov01 <- aov(Zn_liver_mg_kg ~ station, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in Zn_liver_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(Zn_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$Zn_liver_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Zn <- lm(log_prop_hepatovac_area ~ Zn_liver_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Zn)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(Zn_liver_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Zn <- gam(log_prop_hepatovac_area ~ s(Zn_liver_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_Zn)

# Polynomial
poly2_Zn <- lm(log_prop_hepatovac_area ~ poly(Zn_liver_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_Zn)

poly3_Zn <- lm(log_prop_hepatovac_area ~ poly(Zn_liver_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Zn)

### total_metals_mg_kg ------------------------------------------------------
data_clean <- subset(data_TOXEM, !is.na(total_metals_mg_kg))
#data_clean <- na.omit(data_TOXEM[, c("total_metals_mg_kg", "log_prop_hepatovac_area", "age_month", "station", "season","sex", "prop_selected_tissue_area")])

ggplot(data_TOXEM) +
  aes(x = season, y = total_metals_mg_kg, color = season, shape = sex) +
  #geom_boxplot() +
  geom_point() +
  geom_jitter() + 
  facet_wrap(~station)

data_clean %>%
  group_by(station, season, sex) %>%
  summarise(mean_total_metals_mg_kg = mean(total_metals_mg_kg))

aov01 <- aov(total_metals_mg_kg ~ season, data = data_TOXEM)
Anova(aov01, type = "III")
TukeyHSD(aov01, type = "III")
# No significant difference in total_metals_mg_kg level between station, season and sex.

# Linear regression
ggplot(data_TOXEM) +
  aes(total_metals_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = lm,  se = F) +  
  facet_wrap(~ station)

cor.test(x = data_TOXEM$total_metals_mg_kg, y = data_TOXEM$log_prop_hepatovac_area, method = "pearson", use = "complete.obs")

lm_Zn <- lm(log_prop_hepatovac_area ~ total_metals_mg_kg + station + poly(age_month, 2), data = data_TOXEM)
summary(lm_Zn)

# Non-linear regression
# GAM
ggplot(data_TOXEM) +
  aes(total_metals_mg_kg, log_prop_hepatovac_area, col = season, group = station, shape = sex) +
  geom_point() +
  geom_smooth(method = "gam", se = F)

gam_Zn <- gam(log_prop_hepatovac_area ~ s(total_metals_mg_kg) + station + poly(age_month, 2), data = data_TOXEM)
summary(gam_Zn)

# Polynomial
poly2_Zn <- lm(log_prop_hepatovac_area ~ poly(total_metals_mg_kg, 2) + station + poly(age_month, 2), data = data_clean)
summary(poly2_Zn)

poly3_Zn <- lm(log_prop_hepatovac_area ~ poly(total_metals_mg_kg, 3) + station * sex + poly(age_month, 2), data = data_clean)
summary(poly3_Zn)


## plot_metals ------------------------------------------------------------
## plot_PCB -------------------------------------------------
ggarrange(Mn_liver, Cu_liver, Hg_liver, Ni_liver, ncol = 2, nrow = 2, labels = c("A","B","C","D"))

# bkmr_data ---------------------------------------------------------------------
# pollutants sampling design was unbalanced -> not the same number of samples per group for each pollutants
# check number of complete value per group for each pollutants

# Select only pollutants variables
# which(colnames(data_TOXEM)=="naphtalene_liver_ng_g")
# which(colnames(data_TOXEM)=="Mg_liver_mg_kg")

setwd("/home/valentin/Desktop/hepatovac/results/bkmr")

pollutant_data <- data_TOXEM

pollutant_data <- pollutant_data %>% relocate(age_month, .after = season)
# pollutant_data <- pollutant_data %>% relocate(weight_g, .after = season)

which(colnames(pollutant_data)=="naphtalene_liver_ng_g")
which(colnames(pollutant_data)=="Zn_liver_mg_kg")

pollutant_data <- pollutant_data[c(1:7, 73:142)]
pollutant_data <- pollutant_data %>% select(-month)
rownames(pollutant_data) <- pollutant_data[,"image"]
pollutant_data <- pollutant_data %>% select(-image)

colSums(is.na(pollutant_data))
rowSums(is.na(pollutant_data))

# Remove row with more than 71 NAs
pollutant_data <- pollutant_data[!rowSums(is.na(pollutant_data)) >= 71, ]

glimpse(pollutant_data)

## total_data --------------------------------------------------------------
total_data <- pollutant_data %>%
  select(c(station, season, sex, age_month, log_prop_hepatovac_area,
           total_HAP_liver_mg_kg, 
           total_PCB_liver_mg_kg, 
           total_pesticides_liver_mg_kg, 
           total_PBDE_liver_mg_kg,
           Ag_liver_mg_kg,
           Al_liver_mg_kg,
           Cd_liver_mg_kg,
           Co_liver_mg_kg,
           Cr_liver_mg_kg,
           Cu_liver_mg_kg,
           Fe_liver_mg_kg,
           Hg_liver_mg_kg,
           Mn_liver_mg_kg,
           Mo_liver_mg_kg,
           Ni_liver_mg_kg, 
           Pb_liver_mg_kg, 
           Sr_liver_mg_kg, 
           Ti_liver_mg_kg, 
           V_liver_mg_kg, 
           Zn_liver_mg_kg))

# count number of NAs per column
sapply(total_data, function(x) sum(is.na(x)))

# impute missing data with mice package
total_data <- mice(total_data, m = 20, method = 'pmm', seed =123)
total_data <- complete(total_data, 1)

total_data <- total_data %>% select_if(~ !any(is.na(.)))

# recode station and season as numeric variables 
total_data$station <- str_replace_all(total_data$station, 'Canche', '0')
total_data$station <- str_replace_all(total_data$station, 'Seine', '1')
total_data$station <- as.numeric(total_data$station)

total_data$season <- str_replace_all(total_data$season, 'Winter', '0')
total_data$season <- str_replace_all(total_data$season, 'Summer', '1')
total_data$season <- as.numeric(total_data$season)

## reduced_pollutant_data --------------------------------------------------
# Interpolation seems not relevant for bkmr analysis. Thus, we will try to select row and columns with enough complete data for running bkmr.

# remove column with 30% of the values as NAs
reduced_pollutant_data <- pollutant_data[, colSums(!is.na(pollutant_data)) >= 14] # 66% of complete data

# impute missing data with mice package
reduced_pollutant_data <- mice(reduced_pollutant_data, m = 20, method = 'pmm', seed =123)
reduced_pollutant_data <- complete(reduced_pollutant_data, 1)

# replace NAs with median value of the column
#reduced_pollutant_data$column[is.na(reduced_pollutant_data$column)] <- median(reduced_pollutant_data$column, na.rm = TRUE)

reduced_pollutant_data <- reduced_pollutant_data %>% select_if(~ !any(is.na(.)))

reduced_pollutant_data$station <- str_replace_all(reduced_pollutant_data$station, 'Canche', '0')
reduced_pollutant_data$station <- str_replace_all(reduced_pollutant_data$station, 'Seine', '1')
reduced_pollutant_data$station <- as.numeric(reduced_pollutant_data$station)

reduced_pollutant_data$season <- str_replace_all(reduced_pollutant_data$season, 'Winter', '0')
reduced_pollutant_data$season <- str_replace_all(reduced_pollutant_data$season, 'Summer', '1')
reduced_pollutant_data$season <- as.numeric(reduced_pollutant_data$season)

reduced_pollutant_data$sex <- str_replace_all(reduced_pollutant_data$sex, 'M', '0')
reduced_pollutant_data$sex <- str_replace_all(reduced_pollutant_data$sex, 'F', '1')
reduced_pollutant_data$sex <- as.numeric(reduced_pollutant_data$sex)

#reduced_pollutant_data <- reduced_pollutant_data %>% dplyr::select(!starts_with("PCB_"))

## reduced_pollutant_data_canche -------------------------------------------
reduced_pollutant_data_canche <- reduced_pollutant_data %>% filter(station == 0)

## reduced_pollutant_data_seine --------------------------------------------
reduced_pollutant_data_seine <- reduced_pollutant_data %>% filter(station == 1)

## results_list ------------------------------------------------------------
results_list <- list()

for (variable in names(pollutant_data)) {
  if (is.numeric(pollutant_data[[variable]])) {
    temp_data <- pollutant_data %>% drop_na(!!sym(variable)) %>% select_if(~ !any(is.na(.)))  
    results_list[[variable]] <- temp_data
  }
}

# recode categorical variables as numerical factor variables in results_list
for (i in seq_along(results_list)) {
  dataset <- results_list[[i]]
  dataset$station <- str_replace_all(dataset$station, c('Canche' = '0', 'Seine' = '1')) %>%
    as.numeric()
  dataset$season <- str_replace_all(dataset$season, c('Winter' = '0', 'Summer' = '1')) %>%
    as.numeric()
  dataset$sex <- str_replace_all(dataset$sex, c('F' = '0', 'M' = '1')) %>%
     as.numeric()
  results_list[[i]] <- dataset
}

results_list$log_prop_hepatovac_area <- NULL
results_list$age_month <- NULL

# remove dataframe from results_list with less than 6 observations
results_list <- results_list[sapply(results_list, function(x) nrow(x) >= 6)]

rm(i, variable, temp_data, dataset)

# Run bkmr on results_list datasets

# bkmr --------------------------------------------------------------------
# https://jenfb.github.io/bkmr/overview.html
# https://bkmr-guide-iab-env-epi-c1e9f1201284eb8158cc30169fbc7e2f9058900a.gricad-pages.univ-grenoble-alpes.fr/
# https://uncsrp.github.io/TAME2/mixtures-analysis-methods-part-2-bayesian-kernel-machine-regression.html
# https://github.com/Danlu233/Fast_BKMR/blob/main/source/BKMR_function.R

rm(list = ls()[! ls() %in% c("data_TOXEM")])

## total_data ---------------------------------------------------
#exposures <- as.matrix(scale(total_data_data[c(1:10)], center = TRUE))
exposures_total_data <- total_data %>%
  dplyr::select(-c(1:5)) %>%
  scale(center = TRUE) %>%
  as.matrix()

outcome_total_data = total_data |>
  dplyr::select(log_prop_hepatovac_area) |>
  as.matrix()

covariates_total_data = total_data |>
  dplyr::select(c(1:4)) |>
  as.matrix()

cor_matrix <- cor(exposures_total_data)

highly_correlated <- which(abs(cor_matrix) > 0.85 & abs(cor_matrix) < 1 | abs(cor_matrix) < -0.85 & abs(cor_matrix) > -1 , arr.ind = TRUE)
highly_correlated_pairs <- unique(t(apply(highly_correlated, 1, sort)))
if (nrow(highly_correlated_pairs) > 0) {
  for (i in 1:nrow(highly_correlated_pairs)) {
    var1 <- rownames(cor_matrix)[highly_correlated_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[highly_correlated_pairs[i, 2]]
    cat("Highly correlated pair:", var1, "and", var2, " with correlation:", cor_matrix[highly_correlated_pairs[i, 1], highly_correlated_pairs[i, 2]], "\n")
  }}

ncol(exposures_total_data)
which(colnames(exposures_total_data)=="Ag_liver_mg_kg")
which(colnames(exposures_total_data)=="Cd_liver_mg_kg")

fitkm_total_data <- kmbayes(y=outcome_total_data, Z=exposures_total_data, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE,
                            groups = c(1,2,3,4,5,6,5,7,8,9,10,11,12,13,14,15,16,17,18))

saveRDS(fitkm_total_data, file = "fitkm_total_data.rds")
fitkm_total_data <- readRDS("fitkm_total_data.rds")

TracePlot(fit = fitkm_total_data, par = "r", comp = 1) # convergence of r2 (importance of the second exposure)
TracePlot(fit = fitkm_total_data, par = "beta") # convergence of regression coefficients for covariates
TracePlot(fit = fitkm_total_data, par = "sigsq.eps")

# Extract PIPs for main effects
pips_total_data <- ExtractPIPs(fitkm_total_data)
#pips_total_data[order(pips_total_data$PIP),]
pips_total_data[order(pips_total_data$groupPIP),]
# PIPs represents the percentage of iterations when a given exposure was selected.

# Risk overall
risks.overall <- OverallRiskSummaries(fit = fitkm_total_data, y = fitkm_total_data$y, Z = fitkm_total_data$Z, X = fitkm_total_data$X, 
                                      qs = seq(0.1, 0.9, by = 0.10), 
                                      q.fixed = 0.1, method = "exact")

ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  geom_hline(aes(yintercept=0), color="red") +
  labs(y = "Effect on outcome", x="Mixture quantiles", title = "total_data") +
  ggtitle("total_data")

# Exposure-outcome function
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm_total_data, method="exact")

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

# Interactions
risks.singvar <- SingVarRiskSummaries(fit = fitkm_total_data, y = fitkm_total_data$y, Z = fitkm_total_data$Z, X = fitkm_total_data$X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")
risks.singvar

ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()

## reduced_pollutant_data ----------------------------------------------------
exposures = reduced_pollutant_data |>
  select(c(6:27)) |>
  scale() |>
  as.matrix()

outcome = reduced_pollutant_data |>
  select(log_prop_hepatovac_area) |>
  as.matrix()

covariates = reduced_pollutant_data |>
  select(c(1:3)) |>
  as.matrix()

# outcome <- reduced_pollutant_data$log_prop_hepatovac_area
# #outcome <- reduced_pollutant_data$EROD_pmol_min_mg_prot
# 
# exposures <- as.matrix(scale(reduced_pollutant_data[c(7:28)], center = TRUE))
# 
# #covariates <- reduced_pollutant_data$station
# covariates <- as.matrix(reduced_pollutant_data[c(1:4)])

cor_matrix <- cor(exposures)

highly_correlated <- which(abs(cor_matrix) > 0.85 & abs(cor_matrix) < 1 | abs(cor_matrix) < -0.85 & abs(cor_matrix) > -1 , arr.ind = TRUE)
highly_correlated_pairs <- unique(t(apply(highly_correlated, 1, sort)))
if (nrow(highly_correlated_pairs) > 0) {
   for (i in 1:nrow(highly_correlated_pairs)) {
     var1 <- rownames(cor_matrix)[highly_correlated_pairs[i, 1]]
     var2 <- colnames(cor_matrix)[highly_correlated_pairs[i, 2]]
     cat("Highly correlated pair:", var1, "and", var2, " with correlation:", cor_matrix[highly_correlated_pairs[i, 1], highly_correlated_pairs[i, 2]], "\n")
   }}

ncol(exposures)

which(colnames(exposures)=="PCB_49")
which(colnames(exposures)=="PCB_149")

which(colnames(exposures)=="PCB_101")
which(colnames(exposures)=="PCB_118")

#fitkm_reduced_pollutant_data <- kmbayes(y=outcome, Z=exposures, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE)
#fitkm_reduced_pollutant_data <- kmbayes(y=outcome, Z=exposures, X=covariates, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE)
fitkm_reduced_pollutant_data <- kmbayes(y=outcome, Z=exposures, X=covariates, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE,
                                        groups = c(1,2,3,4,5,4,6,2,7, 8, 9, 9,10,11,12,13,14,15,16,17,18,19))
                                       #groups = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)

saveRDS(fitkm_reduced_pollutant_data, file = "fitkm_reduced_pollutant_data.rds")
fitkm_reduced_pollutant_data <- readRDS("fitkm_reduced_pollutant_data.rds")

TracePlot(fit = fitkm_reduced_pollutant_data, par = "r", comp = 1) # convergence of r2 (importance of the second exposure)
TracePlot(fit = fitkm_reduced_pollutant_data, par = "beta") # convergence of regression coefficients for covariates
TracePlot(fit = fitkm_reduced_pollutant_data, par = "sigsq.eps")

# Extract PIPs for main effects
pips_reduced_pollutant_data <- ExtractPIPs(fitkm_reduced_pollutant_data)
#pips_reduced_pollutant_data[order(pips_reduced_pollutant_data$PIP),]
pips_reduced_pollutant_data[order(pips_reduced_pollutant_data$groupPIP),]
# PIPs represents the percentage of iterations when a given exposure was selected.

# Risk overall
risks.overall <- OverallRiskSummaries(fit = fitkm_reduced_pollutant_data, y = fitkm_reduced_pollutant_data$y, Z = fitkm_reduced_pollutant_data$Z, X = fitkm_reduced_pollutant_data$X,
                                      qs = seq(0.1, 0.9, by = 0.10), 
                                      q.fixed = 0.1, method = "exact")

ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  geom_hline(aes(yintercept=0), color="red") +
  labs(y = "Effect on outcome", x="Mixture quantiles", title = "reduced_pollutant_data")

# Exposure-outcome function
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm_reduced_pollutant_data, method="exact")

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

# Interactions
risks.singvar <- SingVarRiskSummaries(fit = fitkm_reduced_pollutant_data, y = fitkm_reduced_pollutant_data$y, Z = fitkm_reduced_pollutant_data$Z, X = fitkm_reduced_pollutant_data$X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")
risks.singvar

ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()

## reduced_pollutant_data_canche ----------------------------------------------------
exposures_canche = reduced_pollutant_data_canche |>
  select(c(6:27)) |>
  scale() |>
  as.matrix()

outcome_canche = reduced_pollutant_data_canche |>
  select(log_prop_hepatovac_area) |>
  as.matrix()

covariates_canche = reduced_pollutant_data_canche |>
  select(c(2:4)) |>
  as.matrix()

# outcome <- reduced_pollutant_data_canche$log_prop_hepatovac_area
# 
# exposures <- as.matrix(scale(reduced_pollutant_data_canche[c(6:27)], center = TRUE))
# 
# covariates <- as.matrix(reduced_pollutant_data_canche[c(2,4)])

cor_matrix <- cor(exposures_canche)

highly_correlated <- which(abs(cor_matrix) > 0.85 & abs(cor_matrix) < 1 | abs(cor_matrix) < -0.85 & abs(cor_matrix) > -1 , arr.ind = TRUE)
highly_correlated_pairs <- unique(t(apply(highly_correlated, 1, sort)))
if (nrow(highly_correlated_pairs) > 0) {
  for (i in 1:nrow(highly_correlated_pairs)) {
    var1 <- rownames(cor_matrix)[highly_correlated_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[highly_correlated_pairs[i, 2]]
    cat("Highly correlated pair:", var1, "and", var2, " with correlation:", cor_matrix[highly_correlated_pairs[i, 1], highly_correlated_pairs[i, 2]], "\n")
  }}

ncol(exposures_canche)

which(colnames(exposures_canche)=="PCB_101")
which(colnames(exposures_canche)=="PCB_118")

#fitkm_reduced_pollutant_data_canche <- kmbayes(y=outcome_canche, Z=exposures_canche, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE)
#fitkm_reduced_pollutant_data_canche <- kmbayes(y=outcome_canche, Z=exposures_canche, X=covariates_canche, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE)
fitkm_reduced_pollutant_data_canche <- kmbayes(y=outcome_canche, Z=exposures_canche, X=covariates_canche, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE,
                                        groups = c(1,2,3,4,5,4,6,7,8, 9,10,11,12,13,14,15,16,17,18,19,20,21))
                                       #groups = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)

saveRDS(fitkm_reduced_pollutant_data_canche, file = "fitkm_reduced_pollutant_data_canche.rds")
#fitkm_reduced_pollutant_data_canche <- readRDS("fitkm_reduced_pollutant_data_canche.rds")

TracePlot(fit = fitkm_reduced_pollutant_data_canche, par = "r", comp = 1) # convergence of r2 (importance of the second exposure)
TracePlot(fit = fitkm_reduced_pollutant_data_canche, par = "beta") # convergence of regression coefficients for covariates
TracePlot(fit = fitkm_reduced_pollutant_data_canche, par = "sigsq.eps")

# Extract PIPs for main effects
pips_reduced_pollutant_data_canche <- ExtractPIPs(fitkm_reduced_pollutant_data_canche)
#pips_reduced_pollutant_data_canche[order(pips_reduced_pollutant_data_canche$PIP),]
pips_reduced_pollutant_data_canche[order(pips_reduced_pollutant_data_canche$groupPIP),]
# PIPs represents the percentage of iterations when a given exposure was selected.

# Risk overall
risks.overall <- OverallRiskSummaries(fit = fitkm_reduced_pollutant_data_canche, y = fitkm_reduced_pollutant_data_canche$y, Z = fitkm_reduced_pollutant_data_canche$Z, X = fitkm_reduced_pollutant_data_canche$X, 
                                      qs = seq(0.1, 0.9, by = 0.10), 
                                      q.fixed = 0.1, method = "exact")

risks.overall$est <- abs(risks.overall$est)

(plot_overall_risk_canche <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  geom_hline(aes(yintercept=0), color="red") +
  labs(y = "Effect on outcome", x="Mixture quantiles", title = "reduced_pollutant_data_canche"))

# Exposure-outcome function
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm_reduced_pollutant_data_canche, method="exact")

pred.resp.univar$est <- abs(pred.resp.univar)

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

# Interactions
risks.singvar <- SingVarRiskSummaries(fit = fitkm_reduced_pollutant_data_canche, y = fitkm_reduced_pollutant_data_canche$y, Z = fitkm_reduced_pollutant_data_canche$Z, X = fitkm_reduced_pollutant_data_canche$X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")

risks.singvar$est <- abs(risks.singvar$est) 

(plot_interactions_canche <- ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip())

## reduced_pollutant_data_seine ----------------------------------------------------
exposures_seine = reduced_pollutant_data_seine |>
  select(c(6:27)) |>
  scale() |>
  as.matrix()

outcome_seine = reduced_pollutant_data_seine |>
  select(log_prop_hepatovac_area) |>
  as.matrix()

covariates_seine = reduced_pollutant_data_seine |>
  select(c(2:4)) |>
  as.matrix()

# outcome <- reduced_pollutant_data_seine$log_prop_hepatovac_area
# 
# exposures <- as.matrix(scale(reduced_pollutant_data_seine[c(6:27)], center = TRUE))
# 
# covariates <- as.matrix(reduced_pollutant_data_seine[c(2:4)])

cor_matrix <- cor(exposures_seine)

highly_correlated <- which(abs(cor_matrix) > 0.85 & abs(cor_matrix) < 1 | abs(cor_matrix) < -0.85 & abs(cor_matrix) > -1 , arr.ind = TRUE)
highly_correlated_pairs <- unique(t(apply(highly_correlated, 1, sort)))
if (nrow(highly_correlated_pairs) > 0) {
  for (i in 1:nrow(highly_correlated_pairs)) {
    var1 <- rownames(cor_matrix)[highly_correlated_pairs[i, 1]]
    var2 <- colnames(cor_matrix)[highly_correlated_pairs[i, 2]]
    cat("Highly correlated pair:", var1, "and", var2, " with correlation:", cor_matrix[highly_correlated_pairs[i, 1], highly_correlated_pairs[i, 2]], "\n")
  }}

ncol(exposures_seine)

which(colnames(exposures_seine)=="PCB_49")
which(colnames(exposures_seine)=="PCB_149")

which(colnames(exposures_seine)=="PCB_101")
which(colnames(exposures_seine)=="PCB_118")

#fitkm_reduced_pollutant_data_seine <- kmbayes(y=outcome, Z=exposures, iter=30000, family = "gaussian", varsel=FALSE, verbose=TRUE)
#fitkm_reduced_pollutant_data_seine <- kmbayes(y=outcome, Z=exposures, X=covariates, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE)
fitkm_reduced_pollutant_data_seine <- kmbayes(y=outcome_seine, Z=exposures_seine, X=covariates_seine, iter=30000, family = "gaussian", varsel=TRUE, verbose=TRUE,
                                        groups = c(1,2,3,4,5,4,6,2,7, 8, 9,10,11,12,13,14,15,16,17,18,19,20))
                                       #groups = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)

saveRDS(fitkm_reduced_pollutant_data_seine, file = "fitkm_reduced_pollutant_data_seine.rds")
#fitkm_reduced_pollutant_data_seine <- readRDS("fitkm_reduced_pollutant_data_seine.rds")

TracePlot(fit = fitkm_reduced_pollutant_data_seine, par = "r", comp = 1) # convergence of r2 (importance of the second exposure)
TracePlot(fit = fitkm_reduced_pollutant_data_seine, par = "beta") # convergence of regression coefficients for covariates
TracePlot(fit = fitkm_reduced_pollutant_data_seine, par = "sigsq.eps")

# Extract PIPs for main effects
pips_reduced_pollutant_data_seine <- ExtractPIPs(fitkm_reduced_pollutant_data_seine)
#pips_reduced_pollutant_data_seine[order(pips_reduced_pollutant_data_seine$PIP),]
pips_reduced_pollutant_data_seine[order(pips_reduced_pollutant_data_seine$groupPIP),]
# PIPs represents the percentage of iterations when a given exposure was selected.

# Risk overall
risks.overall <- OverallRiskSummaries(fit = fitkm_reduced_pollutant_data_seine, y = fitkm_reduced_pollutant_data_seine$y, Z = fitkm_reduced_pollutant_data_seine$Z, X = fitkm_reduced_pollutant_data_seine$X, 
                                      qs = seq(0.1, 0.9, by = 0.10), 
                                      q.fixed = 0.1, method = "exact")

risks.overall$est <- abs(risks.overall$est)

(plot_overall_risk_seine <- ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  geom_hline(aes(yintercept=0), color="red") +
  labs(y = "Effect on outcome", x="Mixture quantiles", title = "reduced_pollutant_data_seine"))

# Exposure-outcome function
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm_reduced_pollutant_data_seine, method="exact")

pred.resp.univar$est <- abs(pred.resp.univar$est)

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

# Interactions
risks.singvar <- SingVarRiskSummaries(fit = fitkm_reduced_pollutant_data_seine, y = fitkm_reduced_pollutant_data_seine$y, Z = fitkm_reduced_pollutant_data_seine$Z, X = fitkm_reduced_pollutant_data_seine$X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")

risks.singvar$est <- abs(risks.singvar$est)

(plot_interactions_seine <- ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip())

## plot_bkmr -------------------------------------------------
ggarrange(plot_overall_risk_seine, plot_overall_risks_canche, ncol = 2, nrow = 1, labels = c("A","B"))
ggarrange(plot_interactions_seine, plot_interactions_canche, ncol = 2, nrow = 1, labels = c("A","B"))

# GLM physio_data_TOXEM ----------------------------------------------------------
# https://www.guru99.com/r-generalized-linear-model.html
physio_data_TOXEM <- subset(data_TOXEM, select = c("log_prop_hepatovac_area", quant_vars_physio_data_TOXEM))

# standardize numeric columns
physio_data_TOXEM <- physio_data_TOXEM %>%
  mutate_if(is.numeric, funs(as.numeric(scale(.))))

# correlation
corr <- data.frame(lapply(physio_data_TOXEM, as.integer))
ggcorr(corr, method = c("pairwise", "spearman"),nbreaks = 6,hjust = 0.8,label = TRUE,label_size = 3,color = "grey50")

# create test & train data
set.seed(1234)
create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1: total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}
data_train <- create_train_test(physio_data_TOXEM, 0.8, train = TRUE)
data_test <- create_train_test(physio_data_TOXEM, 0.8, train = FALSE)
dim(data_train)

# build the model
logit_interaction <- glm(log_prop_hepatovac_area ~ length_cm * weight_g + liver_weight_g + gonad_weight_g + carcasse_weight_g, family = "gaussian", data = data_train)
summary(logit_interaction)
(vif_model <- car::vif(logit_interaction)) # check for multicollinearity -> VIF>10 indicates severe multicollinearity.

logit_interaction_reduced <- glm(log_prop_hepatovac_area ~ length_cm * liver_weight_g + gonad_weight_g, family = "gaussian", data = data_train)
summary(logit_interaction_reduced)
(vif_model <- car::vif(logit_interaction_reduced)) # check for multicollinearity -> VIF>10 indicates severe multicollinearity.

logit_final <- glm(log_prop_hepatovac_area ~ length_cm + liver_weight_g, family = "gaussian", data = data_train)
summary(logit_final)
car::vif(logit_final)  # check for multicollinearity -> VIF>10 indicates severe multicollinearity

# All the variables in physio_data_TOXEM don't show a significant effect on log_prop_hepatovac_area
# Multicollinearity between variables

# GLM biomarker_data_TOXEM ------------------------------------------------
# https://www.guru99.com/r-generalized-linear-model.html
biomarker_data_TOXEM <- subset(data_TOXEM, select = c("log_prop_hepatovac_area", quant_vars_biomarker_data_TOXEM))
summary(biomarker_data_TOXEM)
biomarker_data_TOXEM <- biomarker_data_TOXEM %>% relocate(log_prop_hepatovac_area, .after = prot_muscle_mg_prot_mL_enzymes)

# standardize numeric columns
biomarker_data_TOXEM <- biomarker_data_TOXEM %>%
  mutate_if(is.numeric, funs(as.numeric(scale(.))))

# correlation
corr <- data.frame(lapply(biomarker_data_TOXEM, as.integer))
ggcorr(corr, method = c("pairwise", "spearman"),nbreaks = 6,hjust = 0.8,label = TRUE,label_size = 3,color = "grey50")

# create test & train data
set.seed(1234)
create_train_test <- function(data, size = 0.8, train = TRUE) {
  n_row = nrow(data)
  total_row = size * n_row
  train_sample <- 1: total_row
  if (train == TRUE) {
    return (data[train_sample, ])
  } else {
    return (data[-train_sample, ])
  }
}
data_train <- create_train_test(biomarker_data_TOXEM, 0.8, train = TRUE)
data_test <- create_train_test(biomarker_data_TOXEM, 0.8, train = FALSE)
dim(data_train)

# build the model
logit_interaction <- glm(log_prop_hepatovac_area ~ s(AchE_µmol_min_mg_prot) + 
                           poly(EROD_pmol_min_mg_prot, 2), family = "gaussian", data = data_train)
#poly(TBARS_nmol_eq_MDA_mg_prot, 3) + 
#poly(prot_carbo_nmol_mg, 3) + 
#poly(CS_liver_IU_mg, 2)

summary(logit_interaction)
(vif_model <- car::vif(logit_interaction)) # check for multicollinearity -> VIF>10 indicates severe multicollinearity.

logit_interaction_reduced <- glm(log_prop_hepatovac_area ~ length_cm * liver_weight_g + gonad_weight_g, family = "gaussian", data = data_train)
summary(logit_interaction_reduced)
(vif_model <- car::vif(logit_interaction_reduced)) # check for multicollinearity -> VIF>10 indicates severe multicollinearity.

logit_final <- glm(log_prop_hepatovac_area ~ length_cm + liver_weight_g, family = "gaussian", data = data_train)
summary(logit_final)
car::vif(logit_final)  # check for multicollinearity -> VIF>10 indicates severe multicollinearity

# All the variables in biomarker_data_TOXEM don't show a significant effect on log_prop_hepatovac_area
# Multicollinearity between variables

