rm(list=ls())
#__________________________________________________

# Load raw data of all tetrapod species --> create df_all
# Load threat and habitat information from iucn summary
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

my_path <- "Z:/THESE/5_Data/Distribution_spatiale"

#__________ df_birds with insular & category info ______________

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

setwd("Z:/THESE/6_Projects/predict_vulnerability/Data")
saveRDS(df_birds, "all_birds_MAJ")

##############################################################

setwd("Z:/THESE/6_Projects/predict_vulnerability")

df_birds <- readRDS("Data/all_birds_MAJ")
colSums(is.na(df_birds))

# ------Load threat------
threat_bp <- read.csv("Data/threats_passerif_IUCN_2020_12.csv")
threat_bnp <- read.csv("Data/threats_non_passerif_IUCN_2020_12.csv")
# error in bind_rows: fail to transform code (factor) into numeric
# --> change factor for character
threat_bp$code <- as.character(threat_bp$code)
threat_bnp$code <- as.character(threat_bnp$code)

threat <- bind_rows(threat_bp, threat_bnp)
rm(threat_bnp, threat_bp)


#------- Find associated ias --------

