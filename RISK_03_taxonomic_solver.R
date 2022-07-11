# all sp with native range


rm(list=ls())

library(tidyverse)
library(rredlist)

# open species list 

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
fold <- list.dirs(main_dir)
clean_fold <- fold[-c(1,5)]

clean_fold

all_sp_list <- vector(mode = "list", length = length(clean_fold))

for(i in 1:length(all_sp_list)){
  fil <- list.files(clean_fold[i])
  sp_file <- fil[grepl("sp_list", fil)]
  all_sp_list[[i]] <- readRDS(paste0(clean_fold[i], "/", sp_file))
}

all_sp <- bind_rows(all_sp_list)

length(unique(all_sp$specieskey))
length(unique(all_sp$original_sciname))

# select unique rows without scientific name

all_sp_toclean <- all_sp %>% select(-scientificname) %>% distinct()

# all species original names that are duplicated
dupl_name <- all_sp_toclean$original_sciname[
  duplicated(all_sp_toclean$original_sciname)]

#liste all files in occurrences_clean for selecting files to bind
occ_files <- list.files(
  list.dirs("Output/Occurrences_clean/"),
  full.names = T)

# create df to store new key names
dupl_new_key <- all_sp_toclean %>%
  mutate(new_key = character(nrow(all_sp_toclean)),
         new_species = character(nrow(all_sp_toclean)),
         new_n_occ = character(nrow(all_sp_toclean))) %>%
  filter(original_sciname %in% dupl_name)

for (i in 1:length(dupl_name)){
  #i = 1
  to_deal <- all_sp_toclean %>%
    filter(original_sciname==dupl_name[i])
  
  spk <- unique(to_deal$specieskey)
  
  final_occ <- data.frame()
  final_key = character()
  
  for (j in spk){
    sp_occ <- readRDS(occ_files[grepl(paste0("spk_", j), occ_files)])
    final_occ <- bind_rows(final_occ, sp_occ)
    final_key <- paste0(final_key, j)
  }
  # keep as new name the one with the more occurrences
  new_name <- to_deal$species[which.max(to_deal$n_occ)]
  new_occ <- as.character(nrow(final_occ))
  
  dupl_new_key <- dupl_new_key %>%
    mutate(new_key = if_else(original_sciname==dupl_name[i], final_key, new_key)) %>%
    mutate(new_species = if_else(original_sciname==dupl_name[i], new_name, new_species)) %>%
    mutate(new_n_occ = if_else(original_sciname==dupl_name[i], new_occ, new_n_occ))
  saveRDS(final_occ, 
          paste0("Output/Occurrences_clean_taxo_ok/RISK_03_all_occ_spk_", final_key))
}

dupl_new_key

# attention hemidactylusmercatoryus, platycephalus, mabouia => tout combiner ?
# oui

keys_hemi <- c(5959930, 5959942, 5959946)
final_occ_hemi <- data.frame()
final_key_hemi = character()

for (j in keys_hemi){
  sp_occ <- readRDS(occ_files[grepl(paste0("spk_", j), occ_files)])
  final_occ_hemi <- bind_rows(final_occ_hemi, sp_occ)
  final_key_hemi <- paste0(final_key_hemi, j)
}

# remove the two other files
f_rm_list <- list.files("Output/Occurrences_clean_taxo_ok/",
                        full.names = T)
f_to_rm <- as.list(f_rm_list[grep("5959942", f_rm_list)])
lapply(f_to_rm, file.remove)

# replace with this one
saveRDS(final_occ_hemi, 
        paste0("Output/Occurrences_clean_taxo_ok/RISK_03_all_occ_spk_", final_key_hemi))
dupl_new_key <- dupl_new_key %>%
  mutate(new_key = if_else(new_species=="Hemidactylus mabouia", final_key_hemi, new_key)) %>%
  mutate(new_n_occ = if_else(new_species=="Hemidactylus mabouia", new_occ, new_n_occ))


# get final table of names and species 
all_sp_clean <- bind_rows(
  all_sp_toclean %>%
    mutate(new_key = as.character(specieskey),
           new_species = species,
           new_n_occ = as.character(n_occ)) %>%
    filter(!(original_sciname %in% dupl_name)),
  dupl_new_key) %>%
  distinct(original_sciname, new_key, new_species, new_n_occ) 

# check unicity species_name x key
uni <- all_sp_clean %>%
  distinct(new_key,new_species)

uni$new_key[duplicated(uni$new_key)]
uni$new_species[duplicated(uni$new_species)]

# change Hemidactylus marbouia for new sp key and occ nb
all_sp_clean$new_species[all_sp_clean$new_key=="5959942"] <- new_name
all_sp_clean$new_n_occ[all_sp_clean$new_key=="5959942"] <- new_occ
all_sp_clean$new_key[all_sp_clean$new_key=="5959942"] <- final_key_hemi
# check unicity species_name x key
uni <- all_sp_clean %>%
  distinct(new_key,new_species)
# 327 ias, OK !!


# check that all occurrences are in new directory

# previous dir
occ_files <- list.files(list.dirs("Output/Occurrences_clean/"),  full.names = T)

for (l in 1:nrow(all_sp_clean)){
  #l=1
  files_new_dir <- list.files("Output/Occurrences_clean_taxo_ok/",
                              full.names = T)
  key = paste0("spk_", all_sp_clean$new_key[l])
  
  if(sum(grepl(key, files_new_dir))==0){
    file_to_move <- readRDS(occ_files[grepl(key, occ_files)])
    saveRDS(file_to_move,
            paste0("Output/Occurrences_clean_taxo_ok/RISK_03_all_occ_", key))
  }
}

# how many unique species ?
length(list.files("Output/Occurrences_clean_taxo_ok/", full.names = T))
length(unique(all_sp_clean$new_key))

# 327 unique IAS with GBIF occurrences

# save species list

saveRDS(all_sp_clean, "Output/Occurrences_clean_taxo_ok/RISK_03_ias_list_with_occ")


##### Write corresp table with IUCN native x IAS #####

# load original species list from 
ias_x_native_bmr <- readRDS("characterize_ias_yan/Data/00_IAS_x_native_bmr")

# correct species query
correct_query_gbif <- c(
  "macaca fascicularis ssp. fascicularis" = "Macaca fascicularis subsp. fascicularis", 
  "morelia spilota ssp. imbricata" = "Morelia spilota subsp. imbricata",
  "apis mellifera ssp. scutellata" = "Apis mellifera subsp. scutellata",
  "canis lupus ssp. dingo" ="Canis lupus subsp. dingo",
  "gallus gallus ssp. domesticus" = "Gallus gallus f. domesticus",
  "ovis gmelinii" = "Ovis gmelinii subsp. musimon",
  "kittacincla malabaricus" = "Kittacincla malabarica",
  "anolis wattsi" = "Anolis wattsii")

ias_lower <- unique(ias_x_native_bmr$ias_lower)
sp_list_for_gbif_ok <- ias_lower
sp_list_for_gbif_ok[sp_list_for_gbif_ok %in% names(correct_query_gbif)] <- correct_query_gbif

corresp <- data.frame(ias_lower = ias_lower,
                      original_sciname = sp_list_for_gbif_ok)


all_sp_clean_IUCN <- left_join(
  all_sp_clean,
  corresp, 
  by = "original_sciname")


# bind with ias_x_native
# remove duplicates for the same new species 
# remove rows for ias not in GBIF (because no occurrence, or parasite/bacteria/disease)

ias_x_native_bmr_ias_ok <- left_join(
  ias_x_native_bmr, 
  all_sp_clean_IUCN, 
  by = "ias_lower") %>%
  # select(-c(ias, ias_lower, original_sciname)) %>%
  rename(ias_gbif_key = new_key,
         ias_species = new_species, 
         ias_n_occ_gbif = new_n_occ) %>%
  filter(!is.na(ias_gbif_key)) %>%
  distinct()

  
length(unique(ias_x_native_bmr_ias_ok$scientificName))


saveRDS(ias_x_native_bmr_ias_ok, "Data/RISK_03_ias_x_native_bmr_taxo_ok")

# add taxonomy information

  
gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")



