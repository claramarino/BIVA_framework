# clean species to model
# 140 IAS
# associated native tetrapods


# Focus on species associated to native amniotes
# Count IAS associated to each group/all groups

rm(list=ls())

library(tidyr)
library(dplyr)
library(rredlist)
library(tidyverse)
library(giscoR)
library(rgbif)
library(taxize)

# load ias x native list, from IUCN threat classification scheme
ias_x_native <- readRDS("Output/00_IAS_x_native_species") %>%
  # remove native amphibians
  filter(Class!="Amphibia") %>%
  # remove unspecified ias
  filter(!grepl("unspecified",ias_lower)) %>%
  # remove _old and _new extensions
  mutate(ias_lower = gsub("_old", "", ias_lower)) %>%
  mutate(ias_lower = gsub("_new", "", ias_lower))
sort(unique(ias_x_native$ias_lower))
length(unique(ias_x_native$binomial_iucn))
table(ias_x_native %>% distinct(binomial_iucn, Class) %>% pull (Class))

# load 140 ias to model (more than 1 associated native)
ias145 <- read.csv2("Output/Data_clean/01_145_IAS_to_model.csv")
ias140_gbifkey <- ias145 %>%
  distinct(specieskey, species) %>%
  rename(speciesKey = specieskey,
         taxon = species)
# on tombe ici à 140 espèces uniques, parce que parmi nos 145 il y a des sous espèces
# il faut maintenant retrouver les noms utilisés pour ces espèces dans l'iucn

ias_corresp <- ias145 %>% distinct(species, original_sciname) %>%
  mutate(iucn_name = tolower(original_sciname))
# add species with modified names

# load original 159 spe with 2 or more interactions
names159 <- readRDS("Output/01_names_2_or_more_clean")
names159[!(names159 %in% ias_corresp$iucn_name)]
ias_corresp$iucn_name[!(ias_corresp$iucn_name %in% names159)]

# correct subspeices names following iucn writing 
ias_corresp$iucn_name[ias_corresp$iucn_name == "senegalia catechu"] <- 
  "acacia catechu"
ias_corresp$iucn_name[ias_corresp$iucn_name == "apis mellifera subsp. scutellata"] <- 
  "apis mellifera ssp. scutellata"
ias_corresp$iucn_name[ias_corresp$iucn_name == "canis lupus subsp. dingo"] <- 
  "canis lupus ssp. dingo"
ias_corresp$iucn_name[ias_corresp$iucn_name == "gallus gallus f. domesticus"] <- 
  "gallus gallus ssp. domesticus"
# all sp are now in the list
ias_corresp$iucn_name[!(ias_corresp$iucn_name %in% names159)]
# add sus domesticus
sus <- c("Sus scrofa","sus domesticus", "sus domesticus")
names(sus) = colnames(ias_corresp)
ias_corresp <- bind_rows(ias_corresp, sus)

# 140 species, with their corresponding name in the IUCN table
length(unique(ias_corresp$species))

table(ias145 %>% distinct(class, species) %>% pull(class))


#----- Find all native species associated to the 140 ias ----

ias140_x_native <- left_join(ias_x_native %>% rename(iucn_name = ias_lower), 
                              ias_corresp,
                              by ="iucn_name") %>% filter(!is.na(species))
length(unique(ias140_x_native$species))

native_assoc <- ias140_x_native %>%
  select(binomial:Class) %>%
  distinct()

table(native_assoc$Class)

# how many native are associated to each IAS
nb_nat <- ias140_x_native %>% distinct(binomial, species) %>%
  group_by(species) %>%
  summarise(nb_nat = n())

summary(nb_nat)

# save nb_nat associated for alien attributes in exposure calculation
saveRDS(nb_nat, "Output/focus_ias_range/01_1_nb_assoc_nat_140ias")



#----- species from Wilfried ----

birdnames <- read.csv2("Data/BirdsNames.csv", sep=",")
mamnames <- read.csv2("Data/NAME_mammals_NatComm2019.csv")
amphnames <- read.csv2("Data/AmphibiansNames.csv", sep = ",") %>%
  mutate(binomial = gsub("_", " ", binomial))

# Correspondance with the 140 ias
# only 3 class so not all ias will be inside obviously

table(ias145 %>% distinct(class, species) %>% pull(class))

ias_mam_wt <- ias140_gbifkey %>% filter(taxon %in% mamnames$DistriName)
# 37 sp de mammals sur les 49
ias_bird_wt <- ias140_gbifkey %>% filter(taxon %in% birdnames$binomial)
# 21 sp de birds sur les 22
ias_amph_wt <- ias140_gbifkey %>% filter(taxon %in% amphnames$binomial)
# 2 amphi /2 sp

# mais on n'a pas les autres classes d'ias modélisées (plantes, insectes, reptiles)


# correspondance with the IAS-A
table(native_assoc$Class, native_assoc$insular_endemic)

ass_nat_mam <- native_assoc %>% filter(binomial %in% mamnames$DistriName)
table(ass_nat_mam$insular_endemic)

ass_nat_bird <- native_assoc %>% filter(binomial %in% birdnames$binomial)
table(ass_nat_bird$insular_endemic)


############## Check if the 140 are in GRIIS ##################
griis_gbif_splist <- read_tsv("Data/GRIIS_GBIF_sp_list_download.csv")

griis_focus <- griis_gbif_splist %>% 
  filter(taxonKey %in% ias140_gbif$speciesKey)

# are species in not in dasco in griis ?
occ_dasco <- readRDS("Output/focus_ias_range/01_ias129_occ_dasco_alien_range")
sp_dasco <- unique(occ_dasco$taxon)

no_d <- setdiff(ias140_gbif$taxon, sp_dasco)

no_d %in% griis_focus$species


plot(subset(occ_dasco, taxon =="Felis catus"))

