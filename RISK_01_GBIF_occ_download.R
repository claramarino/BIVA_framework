# get alien species list from ias_x_native_bmr + synonyms
# download occurrences from gbif

rm(list=ls())

library(tidyverse)
library(rredlist)
library(rgbif)
library(stringr)
library(readr)
library(rgdal)
library(spatialEco)
library(pbapply)
library(CoordinateCleaner)

# load species list
ias_x_native_bmr <- readRDS("characterize_ias_yan/Data/00_IAS_x_native_bmr")


########################### PART 1 ####################################

ias_lower <- unique(ias_x_native_bmr$ias_lower)

#### Get taxonomy from gbif for all species #### 

# raw taxonomy
gbif_taxo <- ias_lower %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %>% # combine all data.frames into one
  filter(status != "DOUBTFUL" & matchtype != "FUZZY")
saveRDS(gbif_taxo,"Output/Taxonomy/RISK_01_taxo_gbif_0")

# some species are not found in gbif
setdiff(ias_lower, unique(gbif_taxo$original_sciname))

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

sp_list_for_gbif_ok <- ias_lower

sp_list_for_gbif_ok[sp_list_for_gbif_ok %in% names(correct_query_gbif)] <- correct_query_gbif

# corrected taxonomy
gbif_taxo_1 <- sp_list_for_gbif_ok %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %>% # combine all data.frames into one
  filter(status != "DOUBTFUL" & matchtype != "FUZZY")

saveRDS(gbif_taxo_1,"Output/Taxonomy/RISK_01_taxo_gbif_1")

# check for bacteria, protista, etc. (remove all sp that needs a host)

table(gbif_taxo_1$kingdom)
# no_anim_plant <- gbif_taxo_1 %>% 
#   filter(kingdom %in% c("Bacteria","Chromista","Fungi"))
# # remove species that need a host => keep only plants and animals
# anim_plant_only <- gbif_taxo_1 %>% 
#   filter(kingdom %in% c("Animalia","Plantae"))
# # table(anim_plant_only$phylum)
# # cryptic_phylum <- anim_plant_only %>%
#   filter(phylum %in% c("Ctenophora", "Mollusca", "Nematoda"))
# # keep ctenophora & mollusca but remove nematoda
# check_duplicates <- gbif_taxo_1 %>% 
#   # filter plants and animals
#   filter(kingdom %in% c("Animalia","Plantae")) %>%
#   # remove nematodes
#   filter(phylum != "Nematoda")
# table(check_duplicates$phylum)
# length(unique(check_duplicates$original_sciname))
# # check if some species are synonyms in different kingdoms
# table(check_duplicates %>% distinct(original_sciname, kingdom, phylum) %>%
#         pull(kingdom)) 
# # at least one sp duplicated
# pbm <- check_duplicates %>% distinct(original_sciname, kingdom, phylum)
# # ammophila arenaria in plants and animals 
# # => refer to a plant in IUCN threat details (Thinornis cucullatus)

gbif_taxo_2  <- gbif_taxo_1 %>% 
  # filter plants and animals
  filter(kingdom %in% c("Animalia","Plantae")) %>%
  # remove nematodes
  filter(phylum != "Nematoda") %>%
  # remove one sp of arthropods that is named as a plant (but which is plant in IUCN threat)
  filter(!(original_sciname == "ammophila arenaria" & phylum == "Arthropoda"))
  
check_bis <- gbif_taxo_2 %>% 
  distinct(original_sciname, kingdom, phylum, order, family) # ok
table(check_bis$phylum)

saveRDS(gbif_taxo_2, "Output/Taxonomy/RISK_01_taxo_gbif_2")



########################### PART 2 #####################################
# Download GBIF data

gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")
table(gbif_taxo_2$phylum)
ias_to_check <- unique(gbif_taxo_2$original_sciname)

# fill in your gbif.org credentials 
user <- "claramarino" # your gbif.org username 
pwd <- "data2021" # your gbif.org password
email <- "claramarino665@gmail.com" # your email 


####                     Arthropods                    ####


gbif_arthro <- gbif_taxo_2 %>% filter(phylum == "Arthropoda")

arthro_download <- occ_download(
  pred_in("taxonKey", unique(gbif_arthro$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(arthro_download)
occ_download_get("0368021-210914110416597", path = "Data/GBIF_occurrences")


####                Ctenophora, Mollusca                    ####


gbif_ctenomol <- gbif_taxo_2 %>% 
  filter(phylum == "Ctenophora" | phylum == "Mollusca")

ctenomol_download <- occ_download(
  pred_in("taxonKey", unique(gbif_ctenomol$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(ctenomol_download)
occ_download_get("0371019-210914110416597", path = "Data/GBIF_occurrences")


####                    Tracheophyta                       ####

gbif_plants <- gbif_taxo_2 %>% filter(phylum == "Tracheophyta")

plants_download <- occ_download(
  pred_in("taxonKey", unique(gbif_plants$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(plants_download)

occ_download_get("0371024-210914110416597", path = "Data/GBIF_occurrences")


###                     Chordata - mammals                    ####

gbif_mam <- gbif_taxo_2 %>% filter(class == "Mammalia")

mam_download <- occ_download(
  pred_in("taxonKey", unique(gbif_mam$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(mam_download)

occ_download_get("0371059-210914110416597", path = "Data/GBIF_occurrences")


####                     Chordata - Reptilia                    ####

gbif_rept <- gbif_taxo_2 %>% filter(class == "Reptilia")

rept_download <- occ_download(
  pred_in("taxonKey", unique(gbif_rept$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(rept_download)

occ_download_get("0371070-210914110416597", path = "Data/GBIF_occurrences")


####                     Chordata - Amphibia                    ####

gbif_amph <- gbif_taxo_2 %>% filter(class == "Amphibia")

amph_download <- occ_download(
  pred_in("taxonKey", unique(gbif_amph$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(amph_download)
occ_download_get("0371082-210914110416597", path = "Data/GBIF_occurrences")


####                 Chordata - Actinopterygii                  ####

gbif_fish <- gbif_taxo_2 %>% filter(class == "Actinopterygii")

fish_download <- occ_download(
  pred_in("taxonKey", unique(gbif_fish$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(fish_download)
occ_download_get("0371092-210914110416597", path = "Data/GBIF_occurrences")


# save download metadata

dl_metadata <- list(
  arthro = occ_download_meta("0368021-210914110416597"),
  plants = occ_download_meta(plants_download),
  ctenomol = occ_download_meta(ctenomol_download),
  birds = occ_download_meta(birds_download),
  mam = occ_download_meta(mam_download),
  amph = occ_download_meta(amph_download),
  rept = occ_download_meta(rept_download),
  fish = occ_download_meta(fish_download)
)

saveRDS(dl_metadata, "Data/GBIF_occurrences/metadata_gbif_downloads")


####                     Chordata - birds                      ####

gbif_birds <- gbif_taxo_2 %>% filter(class == "Aves")

# birds have too many occurrences for a single download
# separate in small occurrences number, middle and high

# count number of occ for each sp
df_count <- gbif_birds %>%
  distinct(usagekey, species) %>%
  mutate(nb_occ = numeric(nrow(gbif_birds)))
for (i in 1:nrow(df_count)){
  df_count$nb_occ[i] <- occ_search(taxonKey = df_count$usagekey[i])$meta$count
  print(i)
}

# select sp with a "small" nb of occurrences
low_nb <- df_count %>% filter(nb_occ < 1000000)
sum(low_nb$nb_occ)
# download
low_nb_birds_download <- occ_download(
  pred_in("taxonKey", unique(low_nb$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

meta_low = occ_download_meta(low_nb_birds_download)
occ_download_get(meta_low$key, path = "Data/GBIF_occurrences/Birds/")

# sp with medium nb of occ
mid_nb <- df_count %>% filter(nb_occ > 1000000 & nb_occ< 2500000)
sum(mid_nb$nb_occ)
# download
mid_nb_birds_download <- occ_download(
  pred_in("taxonKey", unique(mid_nb$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

meta_mid = occ_download_meta(mid_nb_birds_download)
occ_download_get(meta_mid$key, path = "Data/GBIF_occurrences/Birds/")


# sp with high nb of occ
high_nb <- df_count %>% filter(nb_occ > 2500000)

# do it one by one?
download <- as.list(high_nb$usagekey)
names(download) = high_nb$usagekey

# Error: A download limitation is exceeded:
#   User claramarino has too many simultaneous downloads; the limit is 3.
# Please wait for some to complete, or cancel any unwanted downloads.  See your user page.

# need to separate by 3 downloads

for (sp in high_nb$usagekey[1:3]){
  occ_download(
    pred_in("taxonKey", sp),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email)
}

list_down <- occ_download_list(user=user,pwd=pwd, limit = 3)
keys <- as.list(list_down$results$key)

meta_high_1 <- lapply(keys, occ_download_meta)
names(meta_high_1) <- high_nb$usagekey[1:3]

for(sp in names(meta_high_1)){
  if(isTRUE(meta_high_1[[sp]]$status == "SUCCEEDED")){
    occ_download_get(meta_high_1[[sp]]$key,
                     path = paste0("Data/GBIF_occurrences/Birds/"))
  }
}

for (sp in high_nb$usagekey[4:6]){
  occ_download(
    pred_in("taxonKey", sp),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user,pwd=pwd,email=email)
}

list_down <- occ_download_list(user=user,pwd=pwd, limit = 3)
keys <- as.list(list_down$results$key)

meta_high_2 <- lapply(keys, occ_download_meta)
names(meta_high_2) <- high_nb$usagekey[4:6]

for(sp in names(meta_high_2)){
  if(isTRUE(meta_high_2[[sp]]$status == "SUCCEEDED")){
    occ_download_get(meta_high_2[[sp]]$key,
                     path = paste0("Data/GBIF_occurrences/Birds/"))
  }
}

# last high nb species
last <- occ_download(
  pred_in("taxonKey", high_nb$usagekey[7]),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)
meta_last <- occ_download_meta(last)
occ_download_get(meta_last$key,
                 path = paste0("Data/GBIF_occurrences/Birds/"))

# save metadata for birds download

dl_metadata_birds <- list(
  birds_low_nb_occ = occ_download_meta(low_nb_birds_download),
  birds_mid_nb_occ = occ_download_meta(mid_nb_birds_download),
  birds_last = meta_last)
dl_metadata_birds <- c(dl_metadata_birds, meta_high_1, meta_high_2)

saveRDS(dl_metadata_birds, "Data/GBIF_occurrences/Birds/metadata_gbif_downloads_birds")



