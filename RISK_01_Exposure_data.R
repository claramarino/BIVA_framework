
rm(list=ls())

library(tidyverse)
library(rredlist)
library(rgbif)
library(stringr)

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


########################### PART 2 ####################################

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



#### Download arthropods occurrences ####


# fill in your gbif.org credentials 
user <- "claramarino" # your gbif.org username 
pwd <- "data2021" # your gbif.org password
email <- "claramarino665@gmail.com" # your email 

arthro_download <- occ_download(
  pred_in("taxonKey", unique(gbif_arthro$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_meta(arthro_download)

occ_download_get("0368021-210914110416597", path = "Data/GBIF_occurrences")
