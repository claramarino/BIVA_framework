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
arthro_occ_clean <- clean_coord_base(arthro_occ)

# remove points with flagged errors
arthro_occ_clean_no_flag <- flag_rm_terr(arthro_occ_clean)

# species difference?
sp_diff(arthro_occ, arthro_occ_clean_no_flag)

final_sp_list_arthro <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(arthro_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Arthropods"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_arthro <- paste0(main_dir, sub_dir1)

final_sp_list_arthro$n_occ <- numeric(nrow(final_sp_list_arthro))
for(spk in unique(arthro_occ_clean_no_flag$speciesKey)){
  occ_spk <- arthro_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_arthro$n_occ[final_sp_list_arthro$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_arthro, "/RISK_02_arthro_occ_spk_", spk))
}

saveRDS(final_sp_list_arthro, paste0(dir_arthro,"/RISK_02_arthro_sp_list_occ_clean"))

#### CTENOMOL ####

ctenomol_key <- dl_metadata$ctenomol$key
# read file
ctenomol_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", ctenomol_key, ".zip")))

# get clean coordinates
ctenomol_occ_clean <- clean_coord_base(ctenomol_occ)

# remove points with flagged errors
ctenomol_occ_clean_no_flag <- flag_rm_terr(ctenomol_occ_clean)

# species difference?
sp_diff(ctenomol_occ, ctenomol_occ_clean_no_flag)

final_sp_list_ctenomol <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(ctenomol_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Cteno_Mol"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_ctenomol <- paste0(main_dir, sub_dir1)

final_sp_list_ctenomol$n_occ <- numeric(nrow(final_sp_list_ctenomol))
for(spk in unique(ctenomol_occ_clean_no_flag$speciesKey)){
  occ_spk <- ctenomol_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_ctenomol$n_occ[final_sp_list_ctenomol$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_ctenomol, "/RISK_02_ctenomol_occ_spk_", spk))
}

saveRDS(final_sp_list_ctenomol, paste0(dir_ctenomol,"/RISK_02_ctenomol_sp_list_occ_clean"))

#### PLANTS ####

plant_key <- dl_metadata$plants$key
# read file
plant_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", plant_key, ".zip")))

# get clean coordinates
plant_occ_clean <- clean_coord_base(plant_occ)

# remove points with flagged errors
plant_occ_clean_no_flag <- flag_rm_terr(plant_occ_clean)

# species difference?
sp_diff(plant_occ, plant_occ_clean_no_flag)

final_sp_list_plant <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(plant_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Plants"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_plant <- paste0(main_dir, sub_dir1)

final_sp_list_plant$n_occ <- numeric(nrow(final_sp_list_plant))
for(spk in unique(plant_occ_clean_no_flag$speciesKey)){
  occ_spk <- plant_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_plant$n_occ[final_sp_list_plant$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_plant, "/RISK_02_plant_occ_spk_", spk))
}

saveRDS(final_sp_list_plant, paste0(dir_plant,"/RISK_02_plant_sp_list_occ_clean"))

length(unique(final_sp_list_plant$specieskey))

#### AMPHIBIANS ####

amph_key <- dl_metadata$amph$key
# read file
amph_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", amph_key, ".zip")))

# get clean coordinates
amph_occ_clean <- clean_coord_base(amph_occ)

# remove points with flagged errors
amph_occ_clean_no_flag <- flag_rm_terr(amph_occ_clean)

# species difference?
sp_diff(amph_occ, amph_occ_clean_no_flag)

final_sp_list_amph <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(amph_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Amphibians"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_amph <- paste0(main_dir, sub_dir1)

final_sp_list_amph$n_occ <- numeric(nrow(final_sp_list_amph))
for(spk in unique(amph_occ_clean_no_flag$speciesKey)){
  occ_spk <- amph_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_amph$n_occ[final_sp_list_amph$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_amph, "/RISK_02_amph_occ_spk_", spk))
}

saveRDS(final_sp_list_amph, paste0(dir_amph,"/RISK_02_amph_sp_list_occ_clean"))


#### REPTILES ####

# get download key = file name
rept_key <- dl_metadata$rept$key
# read file
rept_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", rept_key, ".zip")))

# get clean coordinates
rept_occ_clean <- clean_coord_base(rept_occ)

# remove points with flagged errors
rept_occ_clean_no_flag <- flag_rm_terr(rept_occ_clean)

# species difference?
sp_diff(rept_occ, rept_occ_clean_no_flag)

final_sp_list_rept <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(rept_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Reptiles"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_rept <- paste0(main_dir, sub_dir1)

final_sp_list_rept$n_occ <- numeric(nrow(final_sp_list_rept))
for(spk in unique(rept_occ_clean_no_flag$speciesKey)){
  occ_spk <- rept_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_rept$n_occ[final_sp_list_rept$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_rept, "/RISK_02_rept_occ_spk_", spk))
}

saveRDS(final_sp_list_rept, paste0(dir_rept,"/RISK_02_rept_sp_list_occ_clean"))

#### FISH ####

fish_key <- dl_metadata$fish$key
# read file
fish_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", fish_key, ".zip")))

# get clean coordinates
fish_occ_clean <- clean_coord_base(fish_occ)

# remove points with flagged errors
fish_occ_clean_no_flag <- flag_rm_terr(fish_occ_clean)

# species difference?
sp_diff(fish_occ, fish_occ_clean_no_flag)

final_sp_list_fish <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(fish_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Fish"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_fish <- paste0(main_dir, sub_dir1)

final_sp_list_fish$n_occ <- numeric(nrow(final_sp_list_fish))
for(spk in unique(fish_occ_clean_no_flag$speciesKey)){
  occ_spk <- fish_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_fish$n_occ[final_sp_list_fish$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_fish, "/RISK_02_fish_occ_spk_", spk))
}

saveRDS(final_sp_list_fish, paste0(dir_fish,"/RISK_02_fish_sp_list_occ_clean"))

#### MAMMALS ####

# get download key = file name
mam_key <- dl_metadata$mam$key
# read file
mam_occ <- readr::read_tsv(unzip(paste0("Data/GBIF_occurrences/", mam_key, ".zip")))

# get clean coordinates
mam_occ_clean <- clean_coord_base(mam_occ)

# remove points with flagged errors
mam_occ_clean_no_flag <- flag_rm_terr(mam_occ_clean)

# species difference?
sp_diff(mam_occ, mam_occ_clean_no_flag)

final_sp_list_mam <- gbif_taxo_2 %>%
  filter(specieskey %in% unique(mam_occ_clean_no_flag$speciesKey)) %>%
  dplyr::distinct(specieskey, species, scientificname, original_sciname)

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Mammals"

if (!file.exists(sub_dir1)){
  dir.create(file.path(main_dir, sub_dir1))
}
dir_mam <- paste0(main_dir, sub_dir1)

final_sp_list_mam$n_occ <- numeric(nrow(final_sp_list_mam))
for(spk in unique(mam_occ_clean_no_flag$speciesKey)){
  occ_spk <- mam_occ_clean_no_flag %>%
    filter(speciesKey == spk)
  final_sp_list_mam$n_occ[final_sp_list_mam$specieskey==spk] <- 
    nrow(occ_spk)
  saveRDS(occ_spk, paste0(dir_mam, "/RISK_02_mam_occ_spk_", spk))
}

saveRDS(final_sp_list_mam, paste0(dir_mam,"/RISK_02_mam_sp_list_occ_clean"))

#### BIRDS ####

dl_metadata_b <- readRDS("Data/GBIF_occurrences/Birds/metadata_gbif_downloads_birds")

for (k in 1:length(dl_metadata_b)){
  #k=1

  b_key <- dl_metadata_b[[k]]$key
  
  unzip(paste0("Data/GBIF_occurrences/Birds/", b_key, ".zip"))
  print(paste0("file ", k, " unzipped"))
  
  b_occ <- readr::read_tsv(paste0(b_key, ".csv"))
  print(paste0("file ", k, " open"))
  
  b_occ_clean <- clean_coord_base(b_occ)
  rm(b_occ)
  
  b_occ_clean_no_flag <- flag_rm_terr(b_occ_clean)
  rm(b_occ_clean)
  
  final_sp_list_b <- gbif_taxo_2 %>%
    filter(specieskey %in% unique(b_occ_clean_no_flag$speciesKey)) %>%
    dplyr::distinct(specieskey, species, scientificname, original_sciname)
  
  main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
  sub_dir1 = "Birds"
  
  if (!file.exists(sub_dir1)){
    dir.create(file.path(main_dir, sub_dir1))
  }
  dir_b <- paste0(main_dir, sub_dir1)
  
  final_sp_list_b$n_occ <- numeric(nrow(final_sp_list_b))
  for(spk in unique(b_occ_clean_no_flag$speciesKey)){
    occ_spk <- b_occ_clean_no_flag %>%
      filter(speciesKey == spk)
    final_sp_list_b$n_occ[final_sp_list_b$specieskey==spk] <- 
      nrow(occ_spk)
    saveRDS(occ_spk, paste0(dir_b, "/RISK_02_bird_occ_spk_", spk))
  }
  
  saveRDS(final_sp_list_b, paste0(dir_b,"/temp/RISK_02_bird_sp_list_occ_clean_", k))
  
  
  }


# combine the 9 species list together

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
sub_dir1 = "Birds"

dir_b <- paste0(main_dir, sub_dir1)

sp_list_all <- data.frame()
for (k in 1:9){
  sp_list_k <- readRDS((paste0(dir_b,"/temp/RISK_02_bird_sp_list_occ_clean_", k)))
  sp_list_all <- bind_rows(sp_list_all, sp_list_k)
}


saveRDS(sp_list_all, paste0(dir_b,"/RISK_02_bird_sp_list_occ_clean"))







########################### BROUILLON ####################################


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


# unzip file before
# unzip(paste0("Data/GBIF_occurrences/", birds_key, ".zip"))


# 2 étapes pour les oiseaux
# étape 1 = ouvrir le fichier en chunck, clean coord de base, sauve par espèces 
# étape 2 = appliquer le cleaning avec flag pour chaque sp, 
# puis sauver comme les autres groupes

# étape 1

split = 20
chunk_size <- round(nrow_birds/split, 0)

breaks_start = c((0:(split-1))*chunk_size)

col_names_birds <- colnames(readr::read_tsv(
  paste0(birds_key, ".csv"),
  n_max = 1)) 

for(s in 11:split){
  
  s=11
  library(data.table)
  
  col_names_birds <- fread(paste0(birds_key, ".csv"), nrows = 0, sep = "\t")
  
  b_occ_s <- readr::read_tsv(
    paste0(birds_key, ".csv"),
    col_names = col_names_birds,
    skip = breaks_start[s],
    n_max = chunk_size-1)
  
  print(paste("slot ", s,"/", split, "open"))
  
  # get clean coordinates
  b_occ_clean <- clean_coord_base(b_occ_s)
  rm(b_occ_s)
  
  list_sp_b_occ_clean <- unique(b_occ_clean$speciesKey)
  
  b_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Data/GBIF_occurrences/Birds/"
  
  for(spk in list_sp_b_occ_clean){
    
    occ_sp <- b_occ_clean %>%
      filter(speciesKey == spk) %>%
      mutate(
        taxonKey = as.character(taxonKey),
        speciesKey = as.character(speciesKey),
        coordinatePrecision = as.character(coordinatePrecision),
        year = as.character(year)
      )
    
    if (file.exists(paste0(b_dir,spk))){
      ori_file <- readRDS(paste0(b_dir,spk))
      final_file <- bind_rows(ori_file, occ_sp) %>%
        distinct()
      saveRDS(final_file, paste0(b_dir,spk))
    } else {
      saveRDS(occ_sp, paste0(b_dir,spk))
    }
  }
}


# étape 2







