# RISK find alien range from GRIIS using Anna's script

rm(list=ls())

library(tidyverse)
library(readr)  
library(data.table)


options(stringsAsFactors=FALSE)


# load species list
sp_info <- readRDS("Output/Native_exotic_range/RISK_12_ias_list_with_occ_IUCN_dwld") %>%
  mutate(iucn_output = if_else(is.na(iucn_output), "NA",iucn_output))


## list of alien species to search for
# = all alien species that have no native range from IUCN
SpList <- unique(sp_info$new_species[sp_info$iucn_output == "NA"]) 


## check alien range
path_griis <- "Z:/THESE/5_Data/Alien_data/GRIIS_data_Anna"
list_dir <- list.dirs(path_griis, recursive=FALSE) #list folders' names in GRIIS_data folder - contains all GRIIS datasets except protected areas  


df <- data.frame(matrix(nrow = 0, ncol= 6 ), stringsAsFactors = FALSE) ## create blank dataframe
colnames(df) <- c("scientificName","locationID","countryCode", "occurrenceStatus", "establishmentMeans","Species_or") ## name column of blank df "Species_or" is the original name from the species list  


for (j in 1:length(SpList)){
  
  sp <- SpList[j] 
  
  for (i in 1:length(list_dir)){
    
    if (j == 1){l <- i} else {l <- i+ ((j-1)*221)} #inventory list number so that the rows do not overlap
    
    GRIIS_file_loc <- list_dir[i]
    GRIIS_file <- list.files(paste0(GRIIS_file_loc), recursive=FALSE)
    GRIIS_file <- GRIIS_file[grepl(".csv", GRIIS_file)] #get list of files that end with .csv
    GRIIS_df <- fread(paste0(GRIIS_file_loc,"/",GRIIS_file ),sep = ",") #read file
    
    sp_row <- which(grepl(sp,GRIIS_df$scientificName))
    sp_row1 <- which(grepl(sp,GRIIS_df$acceptedNameUsage))
    scientificName <- which(colnames(GRIIS_df)=="scientificName")
    locationID <- which(colnames(GRIIS_df)=="locationID")
    countryCode <-which(colnames(GRIIS_df)=="countryCode")
    occurrenceStatus <-which(colnames(GRIIS_df)=="occurrenceStatus")
    establishmentMeans <-which(colnames(GRIIS_df)=="establishmentMeans")
    
    if(length(sp_row)>0){
      if(length(sp_row)==1){
        df[l,1] <- GRIIS_df[sp_row,scientificName]
        if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row,locationID]}
        df[l,3] <- GRIIS_df[sp_row,countryCode]
        if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row,occurrenceStatus]}else{df[l,4] <- "NA"}
        df[l,5] <- GRIIS_df[sp_row,establishmentMeans]
        df[l,6] <- sp
      }else
        if(length(sp_row)>1){
          # sp_row2 <- which(GRIIS_df$taxonRank == "SPECIES" & grepl(sp,GRIIS_df$scientificName))
          sp_row2 <- min(sp_row)
          # if(length(sp_row2==0)){sp_row2 <- min(sp_row)}
          df[l,1] <- GRIIS_df[sp_row2,scientificName]
          if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row2,locationID]}
          df[l,3] <- GRIIS_df[sp_row2,countryCode]
          if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row2,occurrenceStatus]}else{df[l,4] <- "NA"}
          df[l,5] <- GRIIS_df[sp_row2,establishmentMeans]
          df[l,6] <- sp
        }
      
    }
    
    if(length(sp_row1)>0){
      if(length(sp_row1)==1){
        df[l,1] <- GRIIS_df[sp_row1,scientificName]
        if(length(locationID)>0){ df[l,2] <- GRIIS_df[sp_row1,locationID]}
        df[l,3] <- GRIIS_df[sp_row1,countryCode]
        if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row1,occurrenceStatus]}else{df[l,4] <- "NA"}
        df[l,5] <- GRIIS_df[sp_row1,establishmentMeans]
        df[l,6] <- sp
      }else
        if(length(sp_row1)>1){
          sp_row3 <- sp_row1[1]
          df[l,1] <- GRIIS_df[sp_row3,scientificName]
          if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row3,locationID]}
          df[l,3] <- GRIIS_df[sp_row3,countryCode]
          if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row3,occurrenceStatus]}else{df[l,4] <- "NA"}
          df[l,5] <- GRIIS_df[sp_row3,establishmentMeans]
          df[l,6] <- sp
        }
    }
    else {next}
    
  }
  
  rm(list=setdiff(ls(), c("df", "SpList", "path_griis"))) ## remove stuff from R memory
  list_dir <- list.dirs(path_griis, recursive=FALSE)
}

df_clean <- na.omit(df) %>% 
  #remove possible duplicates due to scientific name synonyms
  select(-scientificName) %>% distinct()

str(df_clean)

length(unique(df_clean$Species_or))

#### Save country list for each alien species ####

sp_in_griis <- unique(df_clean$Species_or)
# load initial species list for spk
sp_info <- readRDS("Output/Native_exotic_range/RISK_12_ias_list_with_occ_IUCN_dwld")

for(sp in sp_in_griis){
  sp_countries <- df_clean %>% 
    rename(new_species = Species_or) %>%
    filter(new_species == sp)
  sp_key <- unique(sp_info$new_key[sp_info$new_species == sp])
  
  saveRDS(sp_countries, paste0("Output/Native_exotic_range/Exotic_GRIIS",
                               "/RISK_13_exotic_countries_GRIIS_spk_", sp_key))
}


#### Save species list #####

# add specification for species in GRIIS
sp_info_fin <- sp_info %>%
  rename(nat_exo_info = iucn_output) %>%
  mutate(nat_exo_info = if_else(
    new_species %in% sp_in_griis, "EXO_COUNTRY_GRIIS", nat_exo_info))

saveRDS(sp_info_fin, "Output/Native_exotic_range/RISK_13_ias_list_with_occ_IUCN_GRIIS")

table(sp_info_fin %>% distinct(new_key, nat_exo_info) %>% pull(nat_exo_info))




######################## BROUILLON ########################

# phasianus colchicus

data("world")
dat <- world
dat$pc <- ifelse(dat$iso_a2 %in% df_clean$countryCode, "1", "0")
ggplot() +
  geom_sf(data = dat, aes(fill = pc) )

# all occurrences
pc_occ_all <- readRDS("Output/Occurrences_clean_taxo_ok/RISK_03_all_occ_spk_9752149")

pc_occ_all$alien <- ifelse(pc_occ_all$countryCode %in% df_clean$countryCode, "1", "0")
table(pc_occ_all$alien)

library(sf)
pc_alien <- pc_occ_all %>%
  filter(alien=="1") %>%
  mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
  st_as_sf(
    coords = c("LONG", "LAT"),
    agr = "constant",
    crs = "+proj=longlat +datum=WGS84",
    stringsAsFactors = FALSE,
    remove = TRUE)

ggplot() +
  geom_sf(data = dat, aes(fill = pc) ) +
  geom_sf(data = pc_alien)

# castor canadensis

data("world")
dat <- world
dat$cc <- ifelse(dat$iso_a2 %in% df_clean$countryCode, "1", "0")
ggplot() +
  geom_sf(data = dat, aes(fill = pc) )


cc_occ_all <- readRDS("Output/Occurrences_clean_taxo_ok/RISK_03_all_occ_spk_2439838")
cc_occ_all$alien <- ifelse(cc_occ_all$countryCode %in% df_clean$countryCode, "1", "0")
table(cc_occ_all$alien)
pc_alien <- cc_occ_all %>%
  filter(alien=="1") %>%
  mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
  st_as_sf(
    coords = c("LONG", "LAT"),
    agr = "constant",
    crs = "+proj=longlat +datum=WGS84",
    stringsAsFactors = FALSE,
    remove = TRUE)

ggplot() +
  geom_sf(data = dat, aes(fill = cc) ) +
  geom_sf(data = pc_alien)

