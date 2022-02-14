rm(list=ls())
#__________________________________________________

# Load raw data of all tetrapod species --> create df_all
# (using shapefiles from insular tetrapods project)
# Load threat information from IUCN summary
# Deal with synonyms for reptiles sinc db from GARD (not IUCN)

#__________________________________________________

library(tidyr)
library(sp)
library(rgeos)
library(maptools)
library(sf)
library(rgdal)
library(dplyr)
library(rredlist)


##############################################################
# LIST OF ALL NATIVE TETRAPODS FROM IUCN & GARD
##############################################################

my_path <- "Z:/THESE/5_Data/Distribution_spatiale"


#####-----------------------------------------------------
#                       Birds
#####-----------------------------------------------------

# Load raw attribute tables (without shapefiles)

setwd(paste(my_path,"/1_MAJ_BIRDS_Birdlife", sep=""))
birds_all_data <- readRDS("Birds_all_table_MAJ")
birds_contin_data <- readRDS("Birds_contin_all_table_MAJ")

birds_contin <- birds_contin_data %>% 
  distinct(binomial) %>%
  mutate(insular_endemic = 0)

df_birds <- left_join(birds_all_data, birds_contin, by="binomial")
df_birds$insular_endemic[is.na(df_birds$insular_endemic)] <- 1 

df_birds$binomial[df_birds$binomial=="Antigone rubicunda"] <- "Grus rubicunda"
df_birds$binomial[df_birds$binomial=="Antigone canadensis"] <- "Grus canadensis"
df_birds$binomial[df_birds$binomial=="Antigone antigone"] <- "Grus antigone"

# Load the IUCN category 
setwd("Z:/THESE/5_Data/IUCN_summary/MAJ_2020_12/Birds_passerif")
category_pass <- read.csv("simple_summary.csv")
setwd("Z:/THESE/5_Data/IUCN_summary/MAJ_2020_12/Birds_non_passerif")
category_n_pass <- read.csv("simple_summary.csv")

cate_all <- bind_rows(category_pass, category_n_pass)
cate_all <- cate_all %>% 
  select(c(scientificName:genusName, redlistCategory)) %>%
  distinct() %>%
  rename(binomial = scientificName, category = redlistCategory)
str(cate_all)
cate_all$category <- as.factor(cate_all$category)
levels(cate_all$category)
levels(cate_all$category) <- c("CR","DD","EN","EX","EW","LC","NT","VU")

# Add IUCN category to df_birds
df_birds <- left_join(df_birds, cate_all, by="binomial")

# keep presence =1, origin = 1 ou 2, seasonal = 1
df_birds <- df_birds %>% filter(presence==1) %>%
  filter(origin==1 | origin==2) %>%
  filter(seasonal==1)
# remove useless columns and duplicates
df_birds <- df_birds %>% 
  select(binomial, insular_endemic:category) %>%
  mutate(Class="Aves") %>%
  distinct()
table(df_birds$insular_endemic)

setwd("Z:/THESE/6_Projects/predict_vulnerability")
saveRDS(df_birds, "Data/all_birds_MAJ")

#####-----------------------------------------------------
#                       Mammals
#####-----------------------------------------------------

# Continental mammals
setwd(my_path)
mam_contin <- readOGR("0_MAJ_IUCN_2020/Mam_continentaux_2020_12.shp")
mam_contin_data <- mam_contin@data
rm(mam_contin) # remove to save working space

# All mammals
mam_all <- readOGR("0_MAJ_IUCN_2020/Mam_repare.shp")
mam_all_data <- mam_all@data
rm(mam_all)  # remove to save working space

# __________ df_mammals with insular & category info __________
mam_contin <- mam_contin_data %>% 
  distinct(binomial) %>%
  mutate(insular_endemic = 0)

df_mam <- left_join(mam_all_data, mam_contin, by="binomial")
df_mam$insular_endemic[is.na(df_mam$insular_endemic)] <- 1 

# keep presence =1, origin = 1 ou 2, seasonal = 1
df_mam <- df_mam %>% filter(presence==1) %>%
  filter(origin==1 | origin==2) %>%
  filter(seasonal==1)
# retirer les colones inutiles et les doublons
df_mam <- df_mam %>% 
  select(binomial, kingdom:category, insular_endemic) %>%
  mutate(Class="Mammalia") %>%
  distinct()
table(df_mam$insular_endemic)

setwd("Z:/THESE/6_Projects/predict_vulnerability")
saveRDS(df_mam, "Data/all_mammals_MAJ")



#####-----------------------------------------------------
#                      Amphibians
#####-----------------------------------------------------

# Continental amphib
setwd(my_path)
amphi_contin <- readOGR("0_MAJ_IUCN_2020/Amphi_continentaux_2020_12.shp")
amphi_contin_data <- amphi_contin@data
rm(amphi_contin) # remove to save working space

# All amphibians
amphi_all <- readOGR("0_MAJ_IUCN_2020/Amphi_repare.shp")
amphi_all_data <- amphi_all@data
rm(amphi_all)  # remove to save working space

# _________ df_amphi with insular & category info ____________
amphi_contin <- amphi_contin_data %>% 
  distinct(binomial) %>%
  mutate(insular_endemic = 0)

df_amphi <- left_join(amphi_all_data, amphi_contin, by="binomial")
df_amphi$insular_endemic[is.na(df_amphi$insular_endemic)] <- 1 

# keep presence =1, origin = 1 ou 2, seasonal = 1
df_amphi <- df_amphi %>% filter(presence==1) %>%
  filter(origin==1 | origin==2) %>%
  filter(seasonal==1)
# retirer les colones inutiles et les doublons
df_amphi <- df_amphi %>% 
  select(binomial, kingdom:category, insular_endemic) %>%
  mutate(Class="Amphibia") %>%
  distinct()
table(df_amphi$insular_endemic)

setwd("Z:/THESE/6_Projects/predict_vulnerability")
saveRDS(df_amphi, "Data/all_amphibians_MAJ")

#####-----------------------------------------------------
#                       Reptiles
#####-----------------------------------------------------

#   FROM GARD database

# Continental reptiles
setwd(paste(my_path,"/REPTILE_Distrib_GARD.1_dissolved_ranges", sep=""))
rept_contin <- readOGR("Reptiles_GARD_contin.shp")
rept_contin_data <- rept_contin@data
rm(rept_contin) # remove to save working space

# All reptiles
rept_all <- readOGR("Reptiles_GARD_repare.shp")
rept_all_data <- rept_all@data
rm(rept_all)  # remove to save working space

# ________ df_reptiles with insular & category info _________
rept_contin <- rept_contin_data %>% 
  distinct(Binomial) %>%
  mutate(insular_endemic = 0)

# /!\ Make sure data are available for all reptiles
# since the database is not from IUCN --> need to merge with iucn summary?
# no info on seasonality, presence, origin
# see if a lot of reptiles are invasive, using iucn db
# pres_or_sea <- rept_all_data_iucn %>% 
#   distinct(binomial, presence, origin, seasonal) %>%
#   mutate_all(as.factor)
# summary(pres_or_sea)

df_rept <- left_join(rept_all_data, rept_contin, by="Binomial")
df_rept$insular_endemic[is.na(df_rept$insular_endemic)] <- 1 

# no iucn category in GARD
# Load IUCN summary to join with GARD species
setwd("Z:/THESE/5_Data/IUCN_summary/MAJ_2020_12/Amph_Mam_Rept")
category <- read.csv("simple_summary.csv")
category_rept <- category %>% 
  filter(className=="REPTILIA") %>%
  select(scientificName:genusName, redlistCategory) %>%
  distinct() %>%
  rename(Binomial = scientificName, category = redlistCategory)
category_rept$category <- as.factor(category_rept$category)
levels(category_rept$category)
levels(category_rept$category) <- c("CR","DD","EN","EX","EW","LC",
                                    "LC","LC","NT","NT","VU")
# merge with db GARD
df_rept <- left_join(df_rept, category_rept, by="Binomial")
# retirer les colones inutiles et les doublons
df_rept <- df_rept %>% 
  rename(binomial = Binomial) %>%
  select(binomial, insular_endemic:category) %>%
  mutate(Class="Reptilia") %>%
  distinct()
str(df_rept)

colSums(is.na(df_rept))
# pbm of sp without category info --> potentially synonyms?
# use synonyms from insular tetrapods project
# script 00_Extract_attribute_tables_from_shp.R
setwd("Z:/THESE/6_Projects/insular_tetrapods_FD")
# taxize package, itis database 
df_syno_itis <- readRDS("Output/Synonyms/00_synonyms_rept_GARD_itis")
# rredlist package for other synonyms
df_syno_iucn <- readRDS("Output/Synonyms/00_synonyms_rept_GARD_IUCN_MAJ")

# Bind IUCN & ITIS synonyms
# .id & binomial are IUCN binomial names
# for binding with GARD, use synonym column
df_all_syno <- bind_rows(
  df_syno_iucn %>% select(accepted_name, synonym) %>%
    rename(Binomial = accepted_name, binomial = synonym),
  df_syno_itis %>% select(.id, syn_name) %>%
    rename(Binomial = .id, binomial = syn_name)
) %>% distinct()

# Get IUCN category information
syno_cate <- left_join(df_all_syno, category_rept, by="Binomial") %>%
  rename(binomial_iucn = Binomial)

df_rept_no_info <- df_rept %>%
  filter(is.na(category)) %>%
  select(binomial, insular_endemic)

df_rept_suppl <- left_join(df_rept_no_info, syno_cate, by = "binomial") %>%
  filter(!is.na(category))

dupl_bino <- unique(df_rept_suppl[duplicated(df_rept_suppl$binomial),"binomial"])

df_rept_suppl <- df_rept_suppl %>%
  filter(!(binomial %in% dupl_bino)) %>%
  mutate(Class="Reptilia")

# remove sp with no match IUCN, to avoid too much DD species
df_rept <- df_rept %>%
  filter(!(is.na(category)))

# bind df_rept & df with synonyms
df_rept_all <- bind_rows(df_rept, df_rept_suppl) %>%
  mutate(binomial_iucn = if_else(
    is.na(binomial_iucn), binomial, binomial_iucn))

# make sure synonyms are not duplicated
dupl <- df_rept_all[duplicated(df_rept_all$binomial_iucn),"binomial_iucn"]
df_rept_all <- df_rept_all %>%
  filter(!(binomial_iucn %in% dupl))

setwd("Z:/THESE/6_Projects/predict_vulnerability")
saveRDS(df_rept_all, "Data/all_reptiles_GARD_IUCN_MAJ")



##############################################################
# THREATS ASSOCIATED TO NATIVE TETRAPODS FROM IUCN & GARD
##############################################################

rm(list=ls())


setwd("Z:/THESE/6_Projects/predict_vulnerability")

df_amph <- readRDS("Data/all_amphibians_MAJ")
df_birds <- readRDS("Data/all_birds_MAJ")
df_mam <- readRDS("Data/all_mammals_MAJ")
df_rept <- readRDS("Data/all_reptiles_GARD_IUCN_MAJ")

colSums(is.na(df_rept))

# create one table with all tetrapods
df_all <- bind_rows(df_birds %>% 
                      mutate(binomial_iucn = binomial) %>%
                      rename(order = orderName,
                             kingdom = kingdomName,
                             phylum = phylumName,
                             family = familyName, 
                             genus = genusName) %>%
                      select(-className),
                    df_rept%>%
                      rename(order = orderName,
                             kingdom = kingdomName,
                             phylum = phylumName,
                             family = familyName, 
                             genus = genusName) %>%
                      select(-className),
                    df_mam %>% 
                      mutate(binomial_iucn = binomial) %>%
                      rename(order = order_) %>%
                      select(-class),
                    df_amph %>% 
                      mutate(binomial_iucn = binomial) %>%
                      rename(order = order_) %>%
                      select(-class))
6933+10361+5786+7779

# ------Load threat------
threat_amr <- read.csv("Data/threats_amph_mam_rept_IUCN_2020_12.csv") 
threat_bp <- read.csv("Data/threats_passerif_IUCN_2020_12.csv")
threat_bnp <- read.csv("Data/threats_non_passerif_IUCN_2020_12.csv")
# error in bind_rows: fail to transform code (factor) into numeric
# --> change factor for character
threat_amr$code <- as.character(threat_amr$code)
threat_bp$code <- as.character(threat_bp$code)
threat_bnp$code <- as.character(threat_bnp$code)

threat <- bind_rows(threat_amr, threat_bp, threat_bnp)

rm(threat_amr, threat_bnp, threat_bp)

#------- Find associated ias --------

