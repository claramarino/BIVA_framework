# From GBIF data to clean coord
# read gbif tsv files, clean coordinates
# save coord extent = full range for each species
# save species list with coordinates

rm(list=ls())

library(tidyverse)
library(rredlist)
library(rgbif)
library(stringr)
library(readr)
library(spatialEco)
library(CoordinateCleaner)

dl_metadata <- readRDS("Data/GBIF_occurrences/metadata_gbif_downloads")
gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")

# source functions for data cleaning
source("R/clean_coordinates_gbif.R")


#### ARTHROPODS ####

arthro_key <- dl_metadata$arthro$key
# read file
arthro_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", arthro_key, ".zip")))

# get clean coordinates
arthro_occ_clean <- clean_coord(arthro_occ)

# remove points with flagged errors
arthro_occ_clean_no_flag <- flag_rm(arthro_occ_clean)

# species difference?
sp_diff(arthro_occ, arthro_occ_clean_no_flag) #no

# transform to spatial point dataframe
arthro_occ_sp_list <- occ_to_spatial_point(arthro_occ_clean_no_flag)

final_sp_list_arthro <- gbif_taxo_2 %>%
  filter(specieskey %in% names(arthro_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(arthro_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_arthro_full_range")
saveRDS(final_sp_list_arthro, "Output/Occurrences_clean/RISK_01_arthro_sp_list")


#### CTENOMOL ####

# get download key = file name
cm_key <- dl_metadata$ctenomol$key
# read file
cm_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", cm_key, ".zip")))

# get clean coordinates
cm_occ_clean <- clean_coord(cm_occ)

# species difference?
sp_diff(cm_occ, cm_occ_clean) #no

# transform to spatial point dataframe
cm_occ_sp_list <- occ_to_spatial_point(cm_occ_clean)

final_sp_list_cm <- gbif_taxo_2 %>%
  filter(specieskey %in% names(cm_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(cm_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_ctenomol_full_range")
saveRDS(final_sp_list_cm, "Output/Occurrences_clean/RISK_01_ctenomol_sp_list")

#### PLANTS ####

plant_key <- dl_metadata$plants$key
plants_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", plant_key, ".zip")))

# get clean coordinates
plants_occ_clean <- clean_coord(plants_occ)

# species difference?
sp_diff(plants_occ, plants_occ_clean) # 9 species
rm(plants_occ)
# transform to spatial point dataframe
# chunck for large occ files 

chunk_number = 10
obj_to_spl <- unique(plants_occ_clean$speciesKey)
obj_id_spl <- split(obj_to_spl, 
                    cut(seq_along(obj_to_spl), chunk_number, labels = FALSE))

for (chunk in names(obj_id_spl)){
  taxons = obj_id_spl[[chunk]]
  plants_occ_clean_chunk <- plants_occ_clean %>% 
    filter(speciesKey %in% taxons)
  plants_occ_sp_list_chunk <- occ_to_spatial_point(plants_occ_clean_chunk)
  print(chunk)
  saveRDS(plants_occ_sp_list_chunk, 
          paste0("Output/Occurrences_clean/RISK_01_occ_plants_full_range_c", chunk))
}

final_sp_list_plants <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(plants_occ_clean$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(final_sp_list_plants, "Output/Occurrences_clean/RISK_01_plants_sp_list")



#### AMPHIBIANS ####

# get download key = file name
amph_key <- dl_metadata$amph$key
# read file
amph_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", amph_key, ".zip")))

# get clean coordinates
amph_occ_clean <- clean_coord(amph_occ)

# species difference?
sp_diff(amph_occ, amph_occ_clean) #no

# transform to spatial point dataframe
amph_occ_sp_list <- occ_to_spatial_point(amph_occ_clean)

final_sp_list_amph <- gbif_taxo_2 %>%
  filter(specieskey %in% names(amph_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(amph_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_amph_full_range")
saveRDS(final_sp_list_amph, "Output/Occurrences_clean/RISK_01_amph_sp_list")



#### REPTILES ####

# get download key = file name
rept_key <- dl_metadata$rept$key
# read file
rept_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", rept_key, ".zip")))

# get clean coordinates
rept_occ_clean <- clean_coord(rept_occ)

# species difference?
sp_diff(rept_occ, rept_occ_clean) #no

# transform to spatial point dataframe
rept_occ_sp_list <- occ_to_spatial_point(rept_occ_clean)

final_sp_list_rept <- gbif_taxo_2 %>%
  filter(specieskey %in% names(rept_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(rept_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_rept_full_range")
saveRDS(final_sp_list_rept, "Output/Occurrences_clean/RISK_01_rept_sp_list")

#### FISH ####

# get download key = file name
fish_key <- dl_metadata$fish$key
# read file
fish_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", fish_key, ".zip")))

# get clean coordinates
fish_occ_clean <- clean_coord(fish_occ)

# species difference?
sp_diff(fish_occ, fish_occ_clean) #no

# transform to spatial point dataframe
fish_occ_sp_list <- occ_to_spatial_point(fish_occ_clean)

final_sp_list_fish <- gbif_taxo_2 %>%
  filter(specieskey %in% names(fish_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(fish_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_fish_full_range")
saveRDS(final_sp_list_fish, "Output/Occurrences_clean/RISK_01_fish_sp_list")

#### MAMMALS ####

# get download key = file name
mam_key <- dl_metadata$mam$key
# read file
mam_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", mam_key, ".zip")))

# get clean coordinates
mam_occ_clean <- clean_coord(mam_occ)

# species difference?
sp_diff(mam_occ, mam_occ_clean) #no

# transform to spatial point dataframe
mam_occ_sp_list <- occ_to_spatial_point(mam_occ_clean)

final_sp_list_mam <- gbif_taxo_2 %>%
  filter(specieskey %in% names(mam_occ_sp_list)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

saveRDS(mam_occ_sp_list, "Output/Occurrences_clean/RISK_01_occ_mam_full_range")
saveRDS(final_sp_list_mam, "Output/Occurrences_clean/RISK_01_mam_sp_list")

#### BIRDS ####

birds_key <- dl_metadata$birds$key
nrow_birds <- dl_metadata$birds$totalRecords

# unzip file before
# unzip(paste0("Data/GBIF_occurrences/", birds_key, ".zip"))

split = 20
chunk_size <- round(nrow_birds/split, 0)

breaks_start = c((0:(split-1))*chunk_size)

col_names_birds <- colnames(readr::read_tsv(
  paste0(birds_key, ".csv"),
  n_max = 1)) 

for(s in 3:split){
  
  b_occ_s <- readr::read_tsv(
    paste0(birds_key, ".csv"),
    col_names = col_names_birds,
    skip = breaks_start[s],
    n_max = chunk_size-1)
  
  print(paste("slot ", s,"/", split, "open"))
  
  # get clean coordinates
  b_occ_clean <- clean_coord(b_occ_s)
  rm(b_occ_s)
  
  # transform to spatial point dataframe
  b_occ_sp_list <- occ_to_spatial_point(b_occ_clean)
  
  final_sp_list_b <- gbif_taxo_2 %>%
    filter(specieskey %in% names(b_occ_sp_list)) %>%
    dplyr::distinct(specieskey, species, scientificname, original_sciname) %>%
    mutate(chunk = s)
  
  saveRDS(b_occ_sp_list, 
          paste0("Output/Occurrences_clean/RISK_01_occ_birds_full_range_", s))
  saveRDS(final_sp_list_b, 
          paste0("Output/Occurrences_clean/RISK_01_bird_sp_list_", s))
  
}










########################### PART 5 ####################################


gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")
table(gbif_taxo_2$phylum)
ias_to_check <- unique(gbif_taxo_2$original_sciname)

# for each category of IAS, check species data

#-------------- 1. in iucn

# syno_ias <- data.frame()
# for (i in 1:length(ias_lower)){
#   obj <- rl_synonyms(ias_lower[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   syno_ias <- bind_rows(syno_ias, obj$result)
#   saveRDS(syno_ias,"Output/Synonyms/RISK_01_synonyms_ias_356")
# }
iucn_sp <- readRDS("Output/Synonyms/RISK_01_synonyms_ias_356")

#-------------- 2. in gisd

gisd_sp <- read.csv2("Z:/THESE/5_Data/Alien_data/GISD_List_plants_animals_fungi.csv")

#-------------- 3. in First Record Database

frd <- read.csv2("Z:/THESE/5_Data/Alien_data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")

#--------------  3. in checklists from Eduardo

path_checklist <- c("Z:/THESE/5_Data/Alien_data/Species_checklists/")
files = list.files(path = path_checklist)
checklist_eduardo <- lapply(paste0(path_checklist, files), read.csv)
names(checklist_eduardo) <- c("amphi","ants","birds","fish","fungi","mam","plants","rept","spiders")
rm(path_checklist, files)


#### Arthropods ####

gbif_arthro <- gbif_taxo_2 %>% filter(phylum == "Arthropoda")

length(unique(gbif_arthro$original_sciname))

length(unique(gbif_arthro$canonicalname))

gbif_arthro %>% distinct(canonicalname, species, original_sciname)


possib_names <- unique(c(tolower(gbif_arthro$canonicalname), 
                         tolower(gbif_arthro$species), 
                         tolower(gbif_arthro$original_sciname)))

# check if species in iucn
arthro_iucn <- iucn_sp %>% 
  mutate(lower_accepted = tolower(accepted_name),
         lower_syno = tolower(synonym)) %>%
  filter(lower_accepted %in% possib_names |
           lower_syno %in% possib_names)
# aucune

# check if species in GISD
arthro_gisd <- gisd_sp %>%
  mutate(lower_species = tolower(Species)) %>%
  filter(lower_species %in% possib_names)
gisd_tax <- unique(arthro_gisd$lower_species)
# 16 in gisd
no_match_g_arthro <- possib_names[!(possib_names %in% arthro_gisd$lower_species)]
no_match_g_arthro


# check if species in frd

arthro_frd <- frd %>%
  mutate(lower_taxon = tolower(str_replace_all(Taxon, "[^[:alnum:]]", " "))) %>%
  filter(lower_taxon %in% possib_names)

length(unique(arthro_frd$lower_taxon))
frd_tax <- unique(arthro_frd$lower_taxon)
# 20 in FRD
no_match_frd_arthro <- possib_names[!(possib_names %in% arthro_frd$lower_taxon)]
no_match_frd_arthro


setdiff(frd_tax, gisd_tax)
setdiff(gisd_tax, frd_tax)

arthro_frd_or_gisd <- unique(c(frd_tax, gisd_tax))
setdiff(possib_names, arthro_frd_or_gisd)


# check if species in checklist from Eduardo Arlé

spiders = tolower(unique(checklist_eduardo$spiders$gbifDarwinCore))
possib_names %in% spiders # no spiders
ants = tolower(unique(checklist_eduardo$ants$gbifDarwinCore))
possib_names[possib_names %in% ants] # 8 ants



############################ BROUILLON #################################

table(arthro_occ$occurrenceStatus)
table(arthro_occ_clean$coordinatePrecision)
unique(arthro_occ$establishmentMeans)
unique(arthro_occ$establishmentMeans)


# presence only data
# remove specimen from zoo (living_specimen) or museum/collection (preserved_specimen)

sp_before_clean <- unique(arthro_occ$species)

arthro_occ_clean <- arthro_occ %>% 
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION",
                              "OBSERVATION", "OCCURRENCE")) %>%
  filter(occurrenceStatus == "PRESENT")  %>%
  filter(decimalLatitude!= "" | decimalLongitude!= "") %>%
  mutate(LAT = as.numeric(decimalLatitude),
         LONG = as.numeric(decimalLongitude)) %>%
  filter(!(is.na(LAT)|is.na(LONG))) %>%
  dplyr::select(taxonKey, speciesKey, phylum:species, scientificName, 
                LAT, LONG, coordinatePrecision,
                countryCode, year, establishmentMeans)

sp_after_clean <- unique(arthro_occ_clean$species)
setdiff(sp_before_clean, sp_after_clean)
check <- arthro_occ %>% filter(species %in% setdiff(sp_before_clean, sp_after_clean))
# 2 species have record only for preserved specimens

# create a list of spatial point dataframe for each taxonkey (n = 23)
arthro_occ_sp_list <- vector(mode = "list", 
                             length = length(unique(arthro_occ_clean$speciesKey)))
names(arthro_occ_sp_list) <- sort(unique(arthro_occ_clean$speciesKey))

for(sp in names(arthro_occ_sp_list)){
  WGScoor <-  arthro_occ_clean %>%
    filter(speciesKey == sp)
  coordinates(WGScoor)=~LONG+LAT
  proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
  arthro_occ_sp_list[[sp]] <- WGScoor
}



