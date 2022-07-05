# download from gbif occurrences from the 140 IAS

rm(list=ls())

library(rgbif)
library(tidyverse)
library(readr)


# # load 140 ias to model (more than 1 associated native)
# ias145 <- read.csv2("Output/Data_clean/01_145_IAS_to_model.csv")
# ias140_gbifkey <- ias145 %>%
#   distinct(specieskey, species) %>%
#   rename(speciesKey = specieskey,
#          taxon = species)
# # add number off occurrences from gbif
# df_count <- ias140_gbifkey %>% 
#     mutate(nb_occ = numeric(nrow(ias140_gbifkey)))
# 
# for (i in 1:nrow(df_count)){
#     df_count$nb_occ[i] <- occ_search(taxonKey = df_count$speciesKey[i])$meta$count
#     print(i)}
# saveRDS(df_count, "Output/focus_ias_range/01_2_nb_occ_140_ias")

ias140_gbif <- readRDS("Output/focus_ias_range/01_2_nb_occ_140_ias")


# load occurrences in DASCO
occ_dasco <- readRDS("Output/focus_ias_range/01_ias129_occ_dasco_alien_range")


nb_occ_dasco <- occ_dasco@data %>% 
  group_by(taxon) %>% 
  summarise(nb_occ_alien = n())
nb_occ_all <- left_join(ias140_gbif, nb_occ_dasco, by="taxon")

ggplot(nb_occ_all, aes(x = log(nb_occ), y = log(nb_occ_alien))) +
  geom_point()+
  geom_abline(slope = 1, intercept = 0)










