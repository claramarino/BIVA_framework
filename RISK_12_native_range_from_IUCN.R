# Retrieve native range for all species in iucn
# for each sp with occurrences, save a shapefile containing native range


rm(list=ls())

library(tidyverse)
library(stringr)
library(readr)
library(sp)
library(spatialEco)
library(sf)



# load species list
ias_in_iucn <- readRDS("Output/Native_exotic_range/RISK_11_IAS_sp_in_IUCN") %>%
  filter(!is.na(taxonid)) %>% 
  mutate(ias_name = tolower(ias_name)) %>% 
  select(ias_name, scientific_name, kingdom:class, category, 
         population_trend:terrestrial_system) %>%
  distinct()

# load associated key
sp_info <- readRDS("Output/Native_exotic_range/RISK_11_ias_list_with_occ_IUCN_check")


# bind tables
ias_in_iucn_k <- left_join(ias_in_iucn, 
                           sp_info %>% rename(ias_name = spe_lower) %>%
                             distinct(ias_name, new_key), 
                           by = "ias_name")


# load ranges from iucn class by class
# extract species from iucn list 
# extract native range only
# save as unique file per species key
# compile species list found in IUCN at the end

final_sp_list_IUCN <- c()


table(ias_in_iucn$class)

#### AMPHIBIANS ####

# load native ranges from IUCN
amphi_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/0_MAJ_IUCN_2020/Amphi_repare.shp")
# select only the 3 species and native range 
amph_nat_range <- subset(amphi_all, binomial %in% ias_in_iucn$scientific_name &
                    origin %in% c(1,2))
rm(amphi_all)
class(amph_nat_range)

for (sp in unique(amph_nat_range$binomial)){
  sp_poly <- subset(amph_nat_range, binomial == sp)
  sp_key <- ias_in_iucn_k$new_key[ias_in_iucn_k$scientific_name == sp]
  # save native range for eack spk
  saveRDS(sp_poly, paste0("Output/Native_exotic_range/Native_IUCN",
                          "/RISK_12_native_range_IUCN_spk_", sp_key))
}

final_sp_list_IUCN <- c(final_sp_list_IUCN, unique(amph_nat_range$binomial))


#### MAMMALS ####

# load native ranges from IUCN
mam_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/0_MAJ_IUCN_2020/Mam_repare.shp")
# select only the 63 species and native range 
mam_nat_range <- subset(mam_all, binomial %in% ias_in_iucn$scientific_name &
                           origin %in% c(1,2))
rm(mam_all)
head(mam_nat_range)

mam_sp <- unique(mam_nat_range$binomial)
mam_iucn <- ias_in_iucn %>% filter(class=="MAMMALIA")
setdiff(mam_iucn$scientific_name, mam_sp)

# Mus musculus in missing for origin = native => to search with GRIIS?

for (sp in unique(mam_nat_range$binomial)){
  sp_poly <- subset(mam_nat_range, binomial == sp)
  sp_key <- ias_in_iucn_k$new_key[ias_in_iucn_k$scientific_name == sp]
  # save native range for eack spk
  saveRDS(sp_poly, paste0("Output/Native_exotic_range/Native_IUCN",
                          "/RISK_12_native_range_IUCN_spk_", sp_key))
}

length(list.files("Output/Native_exotic_range/Native_IUCN"))

final_sp_list_IUCN <- c(final_sp_list_IUCN, unique(mam_nat_range$binomial))


#### REPTILES ####

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

# is there any missing sp?
b_sp <- unique(bird_nat_range$binomial)
b_iucn <- ias_in_iucn %>% filter(class=="AVES")
setdiff(b_iucn$scientific_name, b_sp) # no !


for (sp in unique(bird_nat_range$binomial)){
  sp_poly <- subset(bird_nat_range, binomial == sp)
  sp_key <- ias_in_iucn_k$new_key[ias_in_iucn_k$scientific_name == sp]
  # save native range for eack spk
  saveRDS(sp_poly, paste0("Output/Native_exotic_range/Native_IUCN",
                          "/RISK_12_native_range_IUCN_spk_", sp_key))
}

final_sp_list_IUCN <- c(final_sp_list_IUCN, unique(bird_nat_range$binomial))

saveRDS(final_sp_list_IUCN, paste0("Output/Native_exotic_range/Native_IUCN",
                                   "/RISK_12_sp_list_native_range_IUCN"))



#### MISSING REPTILES ####

# check in GARD
rept_gard <- st_read("Z:/THESE/5_Data/Distribution_spatiale/REPTILE_Distrib_GARD.1_dissolved_ranges/Reptiles_GARD_repare.shp")

rept_not_in_iucn <- c("Emoia jakati", "Gonatodes albogularis", 
                      "Trachemys scripta", "Phelsuma dubia")
sp_not_iucn <- sp_info %>% filter(ias_in_iucn=="NO") %>% pull(new_species)

rept_not_in_iucn_full <- c(sp_not_iucn, rept_not_in_iucn)


colnames(rept_gard)
rept_gard_nat <- subset(rept_gard, Binomial %in% rept_not_in_iucn_full)

unique(rept_gard_nat$Binomial)

for (sp in unique(rept_gard_nat$Binomial)){
  sp_poly <- subset(rept_gard_nat, Binomial == sp)
  sp_key <- ias_in_iucn_k$new_key[ias_in_iucn_k$scientific_name == sp]
  # save native range for eack spk
  saveRDS(sp_poly, paste0("Output/Native_exotic_range/Native_IUCN",
                          "/RISK_12_native_range_GARD_spk_", sp_key))
}

final_sp_list_IUCN_GARD <- c(final_sp_list_IUCN, unique(rept_gard_nat$Binomial))

saveRDS(final_sp_list_IUCN_GARD, paste0("Output/Native_exotic_range/Native_IUCN",
                                   "/RISK_12_sp_list_native_range_IUCN_GARD"))



