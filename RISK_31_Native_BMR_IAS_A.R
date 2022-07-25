# Retrieve native range for all species in iucn that are IAS-A
# + filter species from ou ias_x_native list


rm(list=ls())

library(tidyverse)
library(stringr)
library(readr)
library(sp)
library(sf)
library(rredlist)

#### Load data and extract ias A/T info ####

# load threats
native_iasa_x_ias <- readRDS("characterize_ias_yan/Data/00_all_native_ias_a")
length(unique(native_iasa_x_ias$scientificName))
table(native_iasa_x_ias %>% distinct(scientificName, className) %>% pull(className))

# compute T, NT, EX-EW col
native_iasa <- native_iasa_x_ias %>%
  distinct(scientificName, className, redlistCategory) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Endangered", "Critically Endangered", "Vulnerable"),
    "T", "NT_DD")) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Extinct", "Extinct in the Wild"), 
    "EX_EW", threatened))


# add a column for species in ias_x_native with specified ias
native_spe_ias <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok") %>%
  distinct(scientificName) %>% pull(scientificName)

native_iasa <- native_iasa %>%
  mutate(ias_modelled = if_else(scientificName %in% native_spe_ias, "YES","NO"))




##### Polygons of IAS A #####

out_fold = "Output/Sensitivity/Polygons_native_IAS_A"

#### MAMMALS

# load native ranges from IUCN
mam_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/0_MAJ_IUCN_2020/Mam_repare.shp")

mam_nat <- mam_all %>% 
  dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% native_iasa$scientificName)%>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry)

length(unique(mam_nat$binomial))
table(native_iasa$className)

native_iasa %>% 
  dplyr::filter(className =="MAMMALIA") %>%
  dplyr::filter(!(scientificName %in% unique(mam_nat$binomial)))
# all missing species are EX_EW

# save native polygons all ias_a mammals
mam_nat_iasa <- left_join(mam_nat, native_iasa %>%
                            select(scientificName, threatened, ias_modelled) %>%
                            rename(binomial = scientificName),
                          by = "binomial")
saveRDS(mam_nat_iasa, paste0(out_fold, "/RISK_31_Polygon_nat_ias_a_MAM"))


#### REPTILES

# load native ranges from IUCN
rept_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/2_MAJ_REPTILES_IUCN_2022/REPTILES.shp")
# select only the 3 species and native range 
rept_nat_range <- subset(rept_all, binomial %in% ias_in_iucn$scientific_name &
                           origin %in% c(1,2))
colnames(rept_all)
class(rept_nat_range)
r_sp <- unique(rept_nat_range$binomial)
r_iucn <- ias_in_iucn %>% filter(class=="REPTILIA")
setdiff(r_iucn$scientific_name, r_sp)

# 4 reptile species are missing
# "Emoia jakati", "Gonatodes albogularis", "Trachemys scripta"    
# "Phelsuma dubia"

for (sp in unique(rept_nat_range$binomial)){
  sp_poly <- subset(rept_nat_range, binomial == sp)
  sp_key <- ias_in_iucn_k$new_key[ias_in_iucn_k$scientific_name == sp]
  # save native range for eack spk
  saveRDS(sp_poly, paste0("Output/Native_exotic_range/Native_IUCN",
                          "/RISK_12_native_range_IUCN_spk_", sp_key))
}

final_sp_list_IUCN <- c(final_sp_list_IUCN, unique(rept_nat_range$binomial))


#### BIRDS ####

# load native ranges from Birdlife
path_birds <- "Z:/THESE/5_Data/Distribution_spatiale/1_MAJ_BIRDS_Birdlife/"
list.files(path_birds)
b1 <- st_read(paste0(path_birds, "Birds_A_2600_repare.shp"))
b1_filt <- subset(b1, binomial %in% ias_in_iucn$scientific_name &
                    origin %in% c(1,2))
rm(b1)
b2 <- st_read(paste0(path_birds, "Birds_A_3400_repare.shp"))
b2_filt <- subset(b2, binomial %in% ias_in_iucn$scientific_name &
                    origin %in% c(1,2))
rm(b2)
b3 <- st_read(paste0(path_birds, "Birds_B_5800_repare.shp"))
b3_filt <- subset(b3, binomial %in% ias_in_iucn$scientific_name &
                    origin %in% c(1,2))
rm(b3)
b4 <- st_read(paste0(path_birds, "Birds_C_5766_repare.shp"))

b4_filt <- subset(b4, binomial %in% ias_in_iucn$scientific_name &
                    origin %in% c(1,2))
rm(b4)

# bind the four files together
bird_nat_range <- bind_rows(b1_filt, b2_filt, b3_filt, b4_filt)



#### Initialize raster grids as in RISK_21 ####
library(raster)

# initialize raster grid
# resolution 0.1
raster0.1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    resolution = 0.1, vals=NULL)
raster0.1[] <- 0
# resolution 1
raster1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                  resolution = 1, vals=NULL)
raster1[] <- 0
