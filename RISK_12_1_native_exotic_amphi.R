# alien range for amphibians


rm(list=ls())

library(tidyverse)
library(stringr)
library(readr)
library(rgdal)
library(spatialEco)
library(rgdal)


# load full occurrences
amph_occ <- readRDS("Output/Occurrences_clean/RISK_01_occ_amph_full_range")
amph_sp <- readRDS("Output/Occurrences_clean/RISK_01_amph_sp_list")

# # load native ranges from IUCN
# amphi_all <- readOGR("Z:/THESE/5_Data/Distribution_spatiale/0_MAJ_IUCN_2020/Amphi_repare.shp")
# 
# # select only the 3 species
# range_3 <- subset(amphi_all, binomial %in% amph_sp$species)
# saveRDS(range_3, "Output/Native_exotic_range/Amphibians/iucn_range")

# exotic occurrences = all occ outside of native range
range_3 <- readRDS("Output/Native_exotic_range/Amphibians/iucn_range")

library(ggplot2)
pt_data = as.data.frame(amph_occ[[1]])
native_fort = fortify(subset(range_3, 
                             origin %in% c(1,2) &
                               binomial == unique(amph_occ[[1]]$species)))
# Plot
ggplot(native_fort, aes(x = long, y = lat, group = group)) +
  geom_point(data = pt_data, aes(x = LONG, y = LAT, group=NA), shape=3, col = "red") + 
  geom_polygon(colour='black', fill='white')



amph_occ_exotic <- amph_occ

spk= names(amph_occ)[[1]]


for (spk in names(amph_occ)){
  species = amph_sp$species[amph_sp$specieskey==spk]
  native_range_sp = subset(range_3, 
                           origin %in% c(1,2) &
                             binomial == species)
  amph_occ_exotic[[spk]] <- erase.point(amph_occ[[spk]], native_range_sp)
}


pt_exo = as.data.frame(amph_occ_exotic[[3]])
pt_native = as.data.frame(amph_occ[[3]])

native_fort = fortify(subset(range_3, 
                             origin %in% c(1,2) &
                               binomial == unique(amph_occ[[3]]$species)))
exotic_fort = fortify(subset(range_3, 
                             origin==3 &
                               binomial == unique(amph_occ[[3]]$species)))
# Plot
ggplot() +
  geom_point(data = pt_native, aes(x = LONG, y = LAT, group=NA), 
             shape=3, col = "red")+
  geom_point(data = pt_exo, aes(x = LONG, y = LAT, group=NA), 
             shape=3, col = "blue") + 
  geom_polygon(data = native_fort, aes(x = long, y = lat, group = group),
               colour='black', fill='white') +
  geom_polygon(data = exotic_fort, aes(x = long, y = lat, group = group),
               colour='black', fill='green')





