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
library(tidyverse)


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
  filter(seasonal %in% c(1:3))
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
  filter(seasonal %in% c(1:3))
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
  filter(seasonal %in% c(1:3))
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

# clean threat data

threat_ias <- threat %>% 
  filter(timing!="Future") %>% # remove future threats
  mutate(code_sub = substr(code, 1, 3)) %>% 
  mutate(code_sub = as.numeric(code_sub),
         all_8 = if_else(code_sub>8 & code_sub<9, 1,0), # check all threat 8
         ias_8.1 = if_else(code_sub==8.1, 1,0)) %>%
  distinct(scientificName, all_8, ias_8.1, ias, stressName, scope, severity, timing)
table(threat_ias$ias_8.1)

ias_to_check <- threat_ias %>% filter(all_8==1 & ias_8.1 ==0)
ias_verif <- threat_ias %>% filter(ias_8.1==1)

# look for species in ias_to_check not present in ias_verif
# permit to see species in cate 8.2/ 8.4 which are not in 8.1
diff <- setdiff(unique(ias_to_check$scientificName) ,
                unique(ias_verif$scientificName))

# select those species to be sure they are not threatened by IAS
ias_to_check_real <- ias_to_check %>% filter(scientificName %in% diff)
# extract all ias names for IAS.8.1 exposed species
threat2 <- threat %>% 
  filter(timing!="Future") %>% # remove future threats
  mutate(code_sub = substr(code, 1, 3)) %>% 
  mutate(code_sub = as.numeric(code_sub),
         all_8 = if_else(code_sub>8 & code_sub<9, 1,0), # check all threat 8
         ias_8.1 = if_else(code_sub==8.1, 1,0))
ias_name <- pull(threat2 %>% filter(ias_8.1==1) %>% distinct(ias), ias)
# remove empty name ""
ias_name <- ias_name[2:length(ias_name)]

ias_to_check <- threat2 %>% filter(all_8==1 & ias_8.1 ==0)
ias_verif <- threat2 %>% filter(ias_8.1==1)
ias_to_check_real <- ias_to_check %>% filter(scientificName %in% diff)

sp_with_named_ias <- unique(pull(ias_to_check_real %>%
                                   filter(ias %in% ias_name), scientificName))

threat_ias <- threat_ias %>% 
  mutate(ias_8.1=if_else(scientificName %in% sp_with_named_ias, 1, ias_8.1))
table(threat_ias$ias_8.1)


#count all ias-t species (but not all have an associated ias)
sp_threat_ias <- threat_ias %>%
  group_by(scientificName) %>%
  summarise(ias_8.1 = sum(ias_8.1)) %>%
  mutate(ias_8.1=if_else(ias_8.1==0,0,1))
table(sp_threat_ias$ias_8.1)
# keep in mind that only tetrapods with details on threat are in this df
# it represents 15 658 sp with IUCN 2020-3 MAJ

#------- Find associated ias --------

# select only 8.1 species (after the cleaning of 8.1 + 8.4 with named ias)

ias_x_native <- left_join(
  threat_ias %>%
    filter(ias_8.1==1) %>%
    filter(ias %in% ias_name) %>%
    rename(binomial_iucn = scientificName),
  df_all, by = "binomial_iucn") %>%
  mutate(
    ias_lower = tolower(ias),
    ias_simple = if_else(grepl("rattus", ias_lower), "rattus spp", ias_lower),
    ias_simple = if_else(grepl("herpest", ias_lower), "herpestidae", ias_simple))

saveRDS(ias_x_native, "Output/00_IAS_x_native_species")

length(unique(ias_x_native$binomial_iucn))
hist(table(ias_x_native$ias))

length(unique(ias_x_native$ias[ias_x_native$Class=="Aves"]))
length(unique(ias_x_native$ias[ias_x_native$Class=="Amphibia"]))
length(unique(ias_x_native$ias[ias_x_native$Class=="Reptilia"]))
length(unique(ias_x_native$ias[ias_x_native$Class=="Mammalia"]))

length(unique(ias_x_native$ias_simple))
dim(ias_x_native %>% distinct(ias_simple, binomial_iucn))

length(unique(ias_x_native$ias_simple[ias_x_native$Class=="Aves"]))
length(unique(ias_x_native$ias_simple[ias_x_native$Class=="Amphibia"]))
length(unique(ias_x_native$ias_simple[ias_x_native$Class=="Reptilia"]))
length(unique(ias_x_native$ias_simple[ias_x_native$Class=="Mammalia"]))


# remove all sp associated to an unspecified ias
# except rattus & herpestes

ias_x_native_spe <- ias_x_native %>%
  filter(!grepl("unspecified",ias_simple))

length(unique(ias_x_native_spe$binomial_iucn))
dim(ias_x_native_spe %>% distinct(ias, binomial_iucn))

length(unique(ias_x_native_spe$ias[ias_x_native_spe$Class=="Aves"]))
length(unique(ias_x_native_spe$ias[ias_x_native_spe$Class=="Amphibia"]))
length(unique(ias_x_native_spe$ias[ias_x_native_spe$Class=="Reptilia"]))
length(unique(ias_x_native_spe$ias[ias_x_native_spe$Class=="Mammalia"]))

length(unique(ias_x_native_spe$ias_simple))
dim(ias_x_native_spe %>% distinct(ias_simple, binomial_iucn))

length(unique(ias_x_native_spe$ias_simple[ias_x_native_spe$Class=="Aves"]))
length(unique(ias_x_native_spe$ias_simple[ias_x_native_spe$Class=="Amphibia"]))
length(unique(ias_x_native_spe$ias_simple[ias_x_native_spe$Class=="Reptilia"]))
length(unique(ias_x_native_spe$ias_simple[ias_x_native_spe$Class=="Mammalia"]))


# Create adjency matrix for each class
# columns = ias, rows = natives
# keep unspecified species for first check

mat_net_list <- vector(mode = "list", length = 4)
names(mat_net_list) <- c("Amphibia","Aves","Mammalia","Reptilia")

# some sp are associated to an ias but don't have a class
# because they are EX/EW and were removed from df_all
dim(ias_x_native_spe %>% filter(is.na(Class)) %>% distinct(binomial_iucn))


for (i in names(mat_net_list)){
  # mat_net_list[[i]] <- ias_x_native %>%
  #   filter(Class == i) %>%
  #   distinct(binomial_iucn, ias_lower) %>%
  #   mutate(count = 1) %>%
  #   pivot_wider(names_from = ias_lower, 
  #               values_from = count, values_fill = 0) %>%
  #   column_to_rownames("binomial_iucn")
  
  # if remove unspecified species

  mat_net_list[[i]] <- ias_x_native_spe %>%
    filter(Class == i) %>%
    distinct(binomial_iucn, ias_simple) %>%
    mutate(count = 1) %>%
    pivot_wider(names_from = ias_simple,
                values_from = count, values_fill = 0) %>%
    column_to_rownames("binomial_iucn")
}

#majority <- c("Majority (50-90%)", "Whole (>90%)")
severe_inter <- c("Very Rapid Declines", "Slow, Significant Declines", 
                  "Rapid Declines", "Causing/Could cause fluctuations")

# significant interactions
ias_x_native_spe_signif <- ias_x_native_spe %>%
  filter(severity %in% severe_inter) #%>%
  #filter(scope %in% majority)
# severity 
# scope

mat_net_list_signif <- vector(mode = "list", length = 4)
names(mat_net_list_signif) <- c("Amphibia","Aves","Mammalia","Reptilia")

for (i in names(mat_net_list_signif)){
  #  remove unspecified species
  # keep only significant interactions
    mat_net_list_signif[[i]] <- ias_x_native_spe_signif %>%
    filter(Class == i) %>%
    distinct(binomial_iucn, ias_simple) %>%
    mutate(count = 1) %>%
    pivot_wider(names_from = ias_simple,
                values_from = count, values_fill = 0) %>%
    column_to_rownames("binomial_iucn")
}

############## Use bipartite package ################

#bipartite
library(bipartite)

adj.matrix <- lapply(mat_net_list_signif, as.matrix)
lapply(adj.matrix, sum) # number of interactions

Sys.time()
# mod_list <- lapply(adj.matrix, computeModules)
# 19 min with scope & severity filtered
# 3h with severity only
Sys.time()
# saveRDS(mod_list, "Output/00_module_bipartite_tetrap_severe")

mod_list <- readRDS("Output/00_module_bipartite_tetrap_severe")

modules_list <- lapply(mod_list, listModuleInformation)

# count number of groups for each class
lapply(modules_list, function(x) length(x[[2]]))

#attribute module to each sp
gps_class <- vector(mode = "list", length = 4)
names(gps_class) <- c("Amphibia","Aves","Mammalia","Reptilia")

for (i in names(mat_net_list_signif)){
  gps_class[[i]] <- ias_x_native_spe_signif %>%
    filter(Class == i) %>%
    filter(category != "EX") %>% # pbm EX EW in reptiles
    filter(category != "EW") %>%
    mutate(module_native = numeric(length(binomial_iucn)),
           module_ias = numeric(length(binomial_iucn)))
  
  for (j in 1:nrow(gps_class[[i]])){
    gps_class[[i]]$module_native[j] = 
      which(grepl(gps_class[[i]]$binomial_iucn[j], modules_list[[i]][[2]]))
    gps_class[[i]]$module_ias[j] = 
      which(grepl(gps_class[[i]]$ias_simple[j], modules_list[[i]][[2]]))
  }
}

# compute network indices for each class
# indices <- lapply(adj.matrix, networklevel) # takes 10 min


# see if groups are reciprocal 
p_class <- lapply(gps_class, function(x){
  ggplot(data = x, aes(x = module_ias, y = module_native)) +
    geom_point(position = "jitter")
})

library(ggpubr)
ggarrange(p_class$Amphibia, p_class$Aves,
          p_class$Mammalia, p_class$Reptilia,
          ncol = 2, nrow = 2)

i="Reptilia" # "Amphibia", "Aves", "Mammalia"
m = as.matrix(table(gps_class[[i]]$module_native, gps_class[[i]]$module_ias))
plotweb(m)


# take into account number of interactions
gps_inter <- lapply(gps_class, function(x){
  nb_inter_ias <- x %>%
    group_by(ias_simple) %>%
    summarize(count_inter_ias = n())
  nb_inter_native <- x %>%
    group_by(binomial_iucn) %>%
    summarize(count_inter_native = n())
  
  gps_inter <- left_join(left_join(x, nb_inter_ias, by="ias_simple"), 
                             nb_inter_native, by = "binomial_iucn")
  return(gps_inter)
})

pias <- lapply(gps_inter, function(x){
  ggplot(data = x, aes(x = module_ias, y = module_native,
                                   color = count_inter_ias)) +
    geom_point(position = "jitter")
})

pnat <- lapply(gps_inter, function(x){
  ggplot(data = x, aes(x = module_ias, y = module_native,
                                   color = count_inter_native)) +
    geom_point(position = "jitter")
})

ggarrange(pias$Amphibia, pnat$Amphibia)
ggarrange(pias$Aves, pnat$Aves)
ggarrange(pias$Mammalia, pnat$Mammalia)
ggarrange(pias$Reptilia, pnat$Reptilia)


# see if less interaction using genus ? families ?

fam_net_list_signif = list()
genus_net_list_signif = list()

for (i in names(mat_net_list_signif)){
  # native genus
  genus_net_list_signif[[i]] <- ias_x_native_spe_signif %>%
    filter(Class == i) %>%
    filter(category != "EX") %>% # pbm EX EW in reptiles
    filter(category != "EW") %>%
    distinct() %>%
    select(genus, ias_simple) %>%
    mutate(count = 1) %>%
    group_by(genus, ias_simple) %>%
    summarise(ninter = sum(count)) %>%
    pivot_wider(names_from = ias_simple,
                values_from = ninter, values_fill = 0) %>%
    column_to_rownames("genus")
  # native family
  fam_net_list_signif[[i]]<- ias_x_native_spe_signif %>%
    filter(Class == i) %>%
    filter(category != "EX") %>% # pbm EX EW in reptiles
    filter(category != "EW") %>%
    distinct() %>%
    select(family, ias_simple) %>%
    mutate(count = 1) %>%
    group_by(family, ias_simple) %>%
    summarise(ninter = sum(count)) %>%
    pivot_wider(names_from = ias_simple,
                values_from = ninter, values_fill = 0) %>%
    column_to_rownames("family")
}



adj.matrix.fam <- lapply(fam_net_list_signif, as.matrix)
lapply(adj.matrix.fam, sum) # number of interactions
plotweb(adj.matrix.fam$Amphibia)
plotweb(adj.matrix.fam$Aves)
plotweb(adj.matrix.fam$Mammalia)
plotweb(adj.matrix.fam$Reptilia)


adj.matrix.gen <- lapply(genus_net_list_signif, as.matrix)
lapply(adj.matrix.gen, sum) # number of interactions



# utilisation du mécanisme

table(gps_inter$Amphibia$stressName)
unique(gps_inter$Aves$stressName)
table(gps_inter$Mammalia$stressName, gps_inter$Mammalia$ias_simple)
table(gps_inter$Reptilia$stressName)

# globalement les IAS sont dans une caté de méca 
# et les méca ne sont pas exclusifs par interaction
# il peut y avoir plusieurs méca pour une interaction native x ias


# est ce que dans les modules les espèces natives sont de même famille ?

table(gps_inter$Amphibia$module_native, gps_inter$Amphibia$family)
table(gps_inter$Aves$module_native, gps_inter$Aves$family)
table(gps_inter$Mammalia$module_native, gps_inter$Mammalia$family)
table(gps_inter$Reptilia$module_native, gps_inter$Reptilia$family)

# est-ce qu'il y a un lien entre module et insular endemic ?
table(gps_inter$Amphibia$module_native, gps_inter$Amphibia$insular_endemic)
table(gps_inter$Aves$module_native, gps_inter$Aves$insular_endemic)
table(gps_inter$Mammalia$module_native, gps_inter$Mammalia$insular_endemic)
table(gps_inter$Reptilia$module_native, gps_inter$Reptilia$insular_endemic)

# test de chi2 (!! wrong test because too small sample for some gps)
chisq.test(table(gps_inter$Amphibia$module_native, gps_inter$Amphibia$insular_endemic))
chisq.test(table(gps_inter$Aves$module_native, gps_inter$Aves$insular_endemic))
chisq.test(table(gps_inter$Mammalia$module_native, gps_inter$Mammalia$insular_endemic))
chisq.test(table(gps_inter$Reptilia$module_native, gps_inter$Reptilia$insular_endemic))


# est-ce qu'il y a un lien entre module et category iucn ?
table(gps_inter$Amphibia$module_native, gps_inter$Amphibia$category)
table(gps_inter$Aves$module_native, gps_inter$Aves$category)
table(gps_inter$Mammalia$module_native, gps_inter$Mammalia$category)
table(gps_inter$Reptilia$module_native, gps_inter$Reptilia$category)

# test avec le module de l'ias
# est-ce qu'il y a un lien entre module ias et insular ?
table(gps_inter$Amphibia$module_ias, gps_inter$Amphibia$insular_endemic)
table(gps_inter$Aves$module_ias, gps_inter$Aves$insular_endemic)
table(gps_inter$Mammalia$module_ias, gps_inter$Mammalia$insular_endemic)
table(gps_inter$Reptilia$module_ias, gps_inter$Reptilia$insular_endemic)


View(gps_inter$Aves)




# nb de natives associées à chaque ias
nb_nat_assoc <- lapply(gps_inter, function(x) {
  x %>% distinct(ias_simple, count_inter_ias)
})
nb_ias_assoc <- lapply(gps_inter, function(x) {
  x %>% distinct(binomial_iucn, count_inter_native)
})

View(nb_ias_assoc$Amphibia)
table(nb_ias_assoc$Amphibia$count_inter_native)
View(nb_ias_assoc$Aves)
table(nb_ias_assoc$Aves$count_inter_native)
View(nb_ias_assoc$Mammalia)
table(nb_ias_assoc$Mammalia$count_inter_native)
View(nb_ias_assoc$Reptilia)
table(nb_ias_assoc$Reptilia$count_inter_native)


View(nb_nat_assoc$Amphibia)
table(nb_nat_assoc$Amphibia$count_inter_ias)
View(nb_nat_assoc$Aves)
table(nb_nat_assoc$Aves$count_inter_ias)
View(nb_nat_assoc$Mammalia)
table(nb_nat_assoc$Mammalia$count_inter_ias)
View(nb_nat_assoc$Reptilia)
table(nb_nat_assoc$Reptilia$count_inter_ias)

############## Use bipartite package ################

#bipartite
library(bipartite)

adj.matrix <- lapply(mat_net_list_signif, as.matrix)
lapply(adj.matrix, sum) # number of interactions

# test with mam only 

g_mam <- as.matrix(mat_net_list$Mammalia)
# mod <- computeModules(g_mam)
# saveRDS(mod, "Output/00_module_bipartite_web_mam")
mod <- readRDS("Output/00_module_bipartite_web_mam")

plotModuleWeb(mod)

module_list = listModuleInformation(mod)
printoutModuleInformation(mod)

View(mod@moduleWeb)
View(mod@originalWeb)
View(mod@modules)
dim(mod@modules)

mod@orderA
mod@orderB
mod@likelihood

module_list[[2]][[1]][[1]]

# attribute to each ias / native its group
gps_mam <- ias_x_native_spe %>%
  filter(Class=="Mammalia") %>%
  mutate(module_native = numeric(nrow(ias_x_native_spe)),
         module_ias = numeric(nrow(ias_x_native_spe)))

for (i in 1:nrow(gps_mam)){
  gps_mam$module_native[i] = which(grepl(gps_mam$binomial_iucn[i], module_list[[2]]))
  gps_mam$module_ias[i] = which(grepl(gps_mam$ias_simple[i], module_list[[2]]))
}
  
# see if groups are reciprocal 
ggplot(data = gps_mam, aes(x = module_ias, y = module_native)) +
  geom_point(position = "jitter")

# take into account number of interactions
nb_inter_ias <- gps_mam %>%
  group_by(ias_simple) %>%
  summarize(count_inter_ias = n())
nb_inter_native <- gps_mam %>%
  group_by(binomial_iucn) %>%
  summarize(count_inter_native = n())

gps_mam_inter <- left_join(left_join(gps_mam, nb_inter_ias, by="ias_simple"), 
                           nb_inter_native, by = "binomial_iucn")

pias <- ggplot(data = gps_mam_inter, aes(x = module_ias, y = module_native,
                                 color = count_inter_ias)) +
  geom_point(position = "jitter")

pnat <- ggplot(data = gps_mam_inter, aes(x = module_ias, y = module_native,
                                 color = count_inter_native)) +
  geom_point(position = "jitter")
library(ggpubr)

ggarrange(pias, pnat)


############## Use blockmodels as in O'connor et al 2020 ################

library(igraph)
library(tidyverse)
library(blockmodels)
library(sna)

my_model <- BM_bernoulli("LBM",M )
my_model$estimate()
which.max(my_model$ICL)

adjm <- adjm.metaweb.diets[, rownames(adjm.metaweb.diets)] ## make sure rows and columns in same order 

adj.matrix <- lapply(mat_net_list, as.matrix) ## convert to matrix format. 

#adj.matrix <- apply(adj.matrix, 2, as.numeric) # make sure elements are numeric!
#rownames(adj.matrix) <- colnames(adj.matrix)  
#diag(adj.matrix) <- 0 ## remove cannibalistic links

lapply(adj.matrix, sum) ## look at number of interactions (if matrix is binary)

# build LBM on adjacency matrix => for bipartite network use "LBM" method
model_list <- lapply(adj.matrix, function(x){
  BM_bernoulli(membership_type = "LBM", adj = x, explore_min=20, explore_max = 30)
})

# explore outputs of model 
# model_list[["Amphibia"]]$estimate()
# model_list[["Aves"]]$estimate()
# model_list[["Mammalia"]]$estimate()
# model_list[["Reptilia"]]$estimate()

lapply(model_list, function(x) x$estimate())

# retrieve best model = which maximizes the ICL 
Qbest <- lapply(model_list, function(x) which.max(x$ICL))
Qbest

plot(model_list[["Amphibia"]]$plot_parameters)

model_list[["Amphibia"]]$memberships[[Qbest[["Amphibia"]]]]$plot()
model_list[["Aves"]]$memberships[[Qbest[["Aves"]]]]$plot()
model_list[["Mammalia"]]$memberships[[Qbest[["Mammalia"]]]]$plot()
model_list[["Reptilia"]]$memberships[[Qbest[["Reptilia"]]]]$plot()


#### species x groups membership table 
spp_groups <- model$memberships[[Qbest]]$Z ## this is the species x groups table (with probabilities of membership)
spp_groups <- as.data.frame(spp_groups)
rownames(spp_groups) <- rownames(adj.matrix)

spp_groups <- Qbest
for (i in 1:length(model_list)){
  spp_groups <- model_list[[i]]$memberships[[Qbest[[i]]]]$Z ## this is the species x groups table (with probabilities of membership)
  spp_groups <- as.data.frame(spp_groups)
  rownames(spp_groups) <- rownames(adj.matrix[[i]])
}

## transform to binary species x group membership
# species belongs to group with maximum probability of membership  (proba ~ 0.99) (else membership = e-05)
colnames(spp_groups) <- 1:Qbest
species.groups <- data.frame(species = rownames(adj.matrix), 
                             group = NA)

rownames(species.groups) <- species.groups$species
for (s in rownames(spp_groups)){
  species.groups[s, "group"] <- as.numeric(which.max(spp_groups[s, ]))
}


###### groups x groups probability interactions = new adjacency matrix at the level of trophic groups 
# Pi = Probability of interactions between classes 
Pi <- model$model_parameters[[Qbest]]$pi

rownames(Pi) <- 1:Qbest
colnames(Pi) <- 1:Qbest

g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
#probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
#distribution
g_mam <- as.matrix(mat_net_list$Mammalia)

#Cluster based on structural equivalence
eq<-equiv.clust(g_mam)
plot(eq)

#Form a blockmodel with distance relaxation of 10
b<-blockmodel(g_mam,eq,h=10, block.content = "density")

b$block.membership
length(unique(b$block.membership))
max(b$block.membership)
# les espèces sont organisées en 35 blocks? 
# attention = dpdt de la distance de relaxation
# en tout il y a 393 identités = 304 natives + 89 IAS

b$block.content
b$order.vector
length(unique(b$order.vector))
304+89

View(b$blocked.data)
View(b$block.model)

b$rlabels
b$cluster.method
b$equiv.fun
b$equiv.metric


plot(b)

