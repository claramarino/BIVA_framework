# exposure completion
# do the current exposure data represent the number of EEE per country?


rm(list=ls())

library(tidyverse)
library(readr)  
library(data.table)
library(sf)
options(stringsAsFactors=FALSE)

# load species list
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")


# path for occurrence files
occ_path <- "Output/True_exposure_alien_species/"
occ_files <- list.files(occ_path)

# initialise empty df
sp_country_55 <- sp_country_110 <- sp_country <-  data.frame()
na_codes <- readRDS(paste0(occ_path, occ_files[29])) %>% filter(is.na(countryCode))

for (i in 1:length(occ_files)){
  print(i)
  
  # load alien occurrences
  occ_sp <- readRDS(paste0(occ_path, occ_files[i]))
  
  # disctinct par country code
  # pour différentes précisions dans les données
  occ_sp_no_filter <- occ_sp %>%
    st_drop_geometry() %>%
    distinct(new_species, new_key, countryCode)
  
  occ_sp_55 <- occ_sp %>% dplyr::filter(coordinateUncertaintyInMeters < 22500) %>%
    st_drop_geometry() %>%
    distinct(new_species, new_key, countryCode)
    
  occ_sp_110 <- occ_sp %>% dplyr::filter(coordinateUncertaintyInMeters < 55000) %>%
    st_drop_geometry() %>%
    distinct(new_species, new_key, countryCode)
  
  # check que toutes les occurrences aient un country code
  if(sum(is.na(occ_sp$countryCode))>0){
    print(paste(i, "contient des NA dans les codes", sep=" "))
    na_codes <- bind_rows(na_codes, 
                          occ_sp %>%
                            filter(is.na(countryCode)))
  }
  
  sp_country <- bind_rows(sp_country, occ_sp_no_filter)
  sp_country_55 <- bind_rows(sp_country_55, occ_sp_55)
  sp_country_110 <- bind_rows(sp_country_110, occ_sp_110)

}

# check for points without country code
world = rnaturalearth::ne_countries(returnclass = "sf")
ggplot(na_codes) +
  geom_sf(data = world)+
  geom_sf()

st_crs(na_codes)
st_crs(world)

# attribute a country code (missing mostly for Namibia)
na_in_country <- st_intersection(na_codes, st_make_valid(world)) 
na_in_country_ok <- na_in_country %>%
  mutate(countryCode = iso_a2) %>% 
  select(new_species:coordinateUncertaintyInMeters) %>%
  st_drop_geometry() %>% distinct()

# add those points in the final list 

sp_country <- bind_rows(
  sp_country,
  na_in_country_ok %>% distinct(new_species, new_key, countryCode)) %>% 
  distinct() %>%
  mutate(InOurDB="yes")
sp_country_55 <- bind_rows(
  sp_country_55,
  na_in_country_ok %>% dplyr::filter(coordinateUncertaintyInMeters < 25500) %>% 
    distinct(new_species, new_key, countryCode)) %>% 
  distinct() %>%
  mutate(InOurDB="yes")
sp_country_110 <- bind_rows(
  sp_country_110,
  na_in_country_ok %>% dplyr::filter(coordinateUncertaintyInMeters < 55000) %>% 
    distinct(new_species, new_key, countryCode)) %>% 
  distinct() %>%
  mutate(InOurDB="yes")


################ species list per country griis ##################


## check alien range
path_griis <- "Z:/THESE/5_Data/Alien_data/GRIIS_data_Anna"
list_dir <- list.dirs(path_griis, recursive=FALSE) #list folders' names in GRIIS_data folder - contains all GRIIS datasets except protected areas  


all_griis <- data.frame()
for (i in 1:length(list_dir)){
  GRIIS_file_loc <- list_dir[i]
  GRIIS_file <- list.files(paste0(GRIIS_file_loc), recursive=FALSE)
  GRIIS_file <- GRIIS_file[grepl(".csv", GRIIS_file)] #get list of files that end with .csv
  GRIIS_df <- fread(paste0(GRIIS_file_loc,"/",GRIIS_file ),sep = ",") %>%
    mutate_all(as.character)
  
  all_griis <- bind_rows(all_griis, GRIIS_df)
  
}

clean_griis <- all_griis %>%
  distinct(scientificName, acceptedNameUsage, locationID, countryCode) %>%
  mutate(scientificName = tolower(scientificName),
         acceptedNameUsage = tolower(acceptedNameUsage))

# take all possible names of ias in our list
ias_sp <- unique(c(sp_info$original_sciname, sp_info$spe_lower))

sp_with_corresp <- data.frame()

for (i in 1:length(ias_sp)){
  match1 <- clean_griis[grep(ias_sp[i], clean_griis$scientificName),]
  match2 <- clean_griis[grep(ias_sp[i], clean_griis$acceptedNameUsage),]
  
  match <- bind_rows(match1, match2) %>%
    mutate(ias_sp_name = ias_sp[i])
  
  sp_with_corresp <- bind_rows(sp_with_corresp, match)
}


length(unique(sp_with_corresp$ias_sp_name))

sp_info_in_griis <- bind_rows(
  left_join(sp_with_corresp, sp_info %>% rename(ias_sp_name = original_sciname)),
  left_join(sp_with_corresp, sp_info %>% rename(ias_sp_name = spe_lower))
) %>% distinct(new_key, new_species, countryCode) %>%
  mutate(InGriis="yes") %>%
  filter(!is.na(new_key)) %>% 
  # some country codes are wrong in GRIIS
  # keep only the two first letters
  mutate(countryCode = strtrim(countryCode, 2)) %>% distinct()


length(unique(sp_info_in_griis$new_key))

# 284 species for which we have a country information

str(sp_country)
sp_country
length(unique(sp_country$new_species))
comp_nolimit <- full_join(sp_country, sp_info_in_griis) %>%
  # keep only the 284 sp for which we have info in GRIIS
  filter(new_species %in% sp_info_in_griis$new_species) %>%
  replace_na(list(InOurDB = "no", InGriis = "no")) %>%
  mutate(cate = case_when(
    InOurDB == "yes" & InGriis =="yes" ~ "In our DB & GRIIS",
    InOurDB == "no" & InGriis =="yes" ~ "Only GRIIS",
    InOurDB == "yes" & InGriis =="no" ~ "Only our DB"
  )) %>% group_by(countryCode, cate) %>%
  summarise(nb_sp = n())

tot <- comp_nolimit %>% group_by(countryCode) %>%
  summarise(tot_sp_count = sum(nb_sp))

prop_country <- left_join(comp_nolimit, tot)%>%
  mutate(prop = nb_sp/tot_sp_count)

table(comp_55$InOurDB, comp_55$InGriis)

comp_55 <- full_join(sp_country_55, sp_info_in_griis) %>%
  replace_na(list(InOurDB = "no", InGriis = "no"))

comp_110 <- full_join(sp_country_110, sp_info_in_griis) %>%
  replace_na(list(InOurDB = "no", InGriis = "no"))





############## Map countries completion ###################

poly_iso <- st_read("Z:/THESE/5_Data/Distribution_spatiale/Countries_ISO3166_boundaries/iso3166-1-boundaries.shp")
str(poly_iso)

ggplot(poly_iso) + geom_sf()

poly_prop <- left_join(poly_iso %>% rename(countryCode = ISO3166.1), 
                      prop_country %>% select(countryCode, tot_sp_count, cate, prop) %>%
                        pivot_wider(names_from = cate, values_from = prop))

colnames(poly_prop)

# prop in griis & our db
ggplot(poly_prop) + geom_sf(aes(fill=`In our DB & GRIIS`), color = NA)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar") +
  geom_sf_text(aes(label=tot_sp_count), color ="white", size =2)

ggplot(poly_prop) + geom_sf(aes(fill=`Only our DB`), color = NA)+
  scale_fill_viridis_c(option = "C")

ggplot(poly_prop) + geom_sf(aes(fill=`Only GRIIS`), color = NA)+
  scale_fill_viridis_c(option = "C") +
  geom_sf_text(aes(label=tot_sp_count))
