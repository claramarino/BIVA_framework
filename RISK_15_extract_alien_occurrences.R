# export only alien points for each sp
# create a df as in DASCO  of alien points for all species in the list 
# species, original name, new_key, x, y, country?


rm(list=ls())

library(tidyverse)
library(stringr)
library(readr)
library(sp)
library(sf)

# species key
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")


# file list occurrence points
occ_path <- "Output/Occurrences_clean_taxo_ok/"
occ_files <- list.files(occ_path)


# file list output
# do only sp for which you haven't done it yet
out_path <- "Output/True_exposure_alien_species/"
# don't forget to clean folder if doing again the procedure
# out_files <- list.files(out_path)
out_files <- character()

##### 1. Native range from IUCN (n = 184 sp) ####

# file list native ranges
nat_path_iucn <- "Output/Native_exotic_range/Native_IUCN/"
nat_files_iucn <- list.files(nat_path_iucn)

sp_key_nat_range_iucn <- sp_info %>% 
  filter(nat_exo_info %in% c("NAT_RANGE_POLY")) %>%
  distinct(new_key) %>% pull(new_key)

# extract exotic occ for each sp in the list
for (k in sp_key_nat_range_iucn){
  # k = "2498305"
  nat_k <- nat_files_iucn[grepl(paste0("spk_", k), nat_files_iucn)]
  out_k <- out_files[grepl(paste0("spk_", k), out_files)]
  print(k)
  
  if(is_empty(nat_k)){ next }
  if(is_empty(out_k)) {
    nat_range <- readRDS(paste0(nat_path_iucn, nat_k))
    
    occ_k <- occ_files[grepl(paste0("spk_", k), occ_files)]
    occ <- readRDS(paste0(occ_path, occ_k))
    
    # convert occ data frame into sf object
    occ_sf <- occ %>%
      mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
      st_as_sf(
        coords = c("LONG", "LAT"),
        agr = "constant",
        crs = "+proj=longlat +datum=WGS84",
        stringsAsFactors = FALSE,
        remove = TRUE)
    
    # create a control column for removing occurrences inside poly
    nat_range <- nat_range %>%
      mutate(control = 1:nrow(nat_range))
    
    # repair potential invalid geometry
    nat_range <- st_make_valid(nat_range)
    precision = 1e7
    
    # if geometry still invalid, decrease precision
    while (sum(st_is_valid(nat_range))<nrow(nat_range)){
      precision = precision/2
      nat_range <- st_make_valid(st_set_precision(nat_range, precision))
    }
    
    # add a buffer around nat range for covering points at the limit
    #plot(nat_range%>%select(control))
    nat_range_buff <- st_buffer(nat_range, dist = 50000)
    #plot(nat_range_buff%>%select(control))
    
    # join points and polygons 
    occ_in_nat_buff <- st_join(occ_sf, nat_range_buff, join = st_within, left = T)

    # select only points outside polygons
    occ_out_nat_buff <- occ_in_nat_buff %>%
      filter(is.na(control)) %>%
      mutate(new_key = k, 
             new_species = unique(sp_info$new_species[sp_info$new_key==k])) %>%
      select(new_species, new_key, countryCode, geometry)
    
    # plot(occ_out_nat %>% select(new_species))
    # plot(occ_out_nat_buff %>% select(new_species))
    
    # save in final data for next analysis
    # some species might have no record outside native range
    # thus save output only if at least one row in occ_out_nat_buff
    
    if(nrow(occ_out_nat_buff)>0){
      saveRDS(occ_out_nat_buff, paste0(out_path, "RISK_15_alien_occ_IUCN_spk_", k))
    }
    print(k)
  }
  
}

length(list.files(out_path))

##### 2. Native range from Other (CABI/Hand) (n = 7 sp) ####

# file list native ranges
nat_path_other <- "Output/Native_exotic_range/Other_database/"
nat_files_other <- list.files(nat_path_other)

sp_key_nat_range_other <- sp_info %>% 
  filter(nat_exo_info %in% c("NAT_COUNTRY_CABI", "NAT_COUNTRY_HAND")) %>%
  distinct(new_key) %>% pull(new_key)

# extract exotic occ for each sp in the list
for (k in sp_key_nat_range_other){
  # k = ""
  nat_k <- nat_files_other[grepl(paste0("spk_", k), nat_files_other)]
  out_k <- out_files[grepl(paste0("spk_", k), out_files)]
  print(k)
  
  if(is_empty(nat_k)){ next }
  if(is_empty(out_k)) {
    nat_range <- readRDS(paste0(nat_path_other, nat_k))
    
    occ_k <- occ_files[grepl(paste0("spk_", k), occ_files)]
    occ <- readRDS(paste0(occ_path, occ_k))
    
    # convert occ data frame into sf object
    occ_sf <- occ %>%
      mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
      st_as_sf(
        coords = c("LONG", "LAT"),
        agr = "constant",
        crs = "+proj=longlat +datum=WGS84",
        stringsAsFactors = FALSE,
        remove = TRUE)
    
    # create a control column for removing occurrences inside poly
    nat_range <- nat_range %>%
      mutate(control = 1:nrow(nat_range))
    
    # repair potential invalid geometry
    nat_range <- st_make_valid(nat_range)
    precision = 1e7
    
    # if geometry still invalid, decrease precision
    while (sum(st_is_valid(nat_range))<nrow(nat_range)){
      precision = precision/2
      nat_range <- st_make_valid(st_set_precision(nat_range, precision))
    }
    
    # join points and polygons 
    occ_in_nat <- st_join(occ_sf, nat_range, join = st_within, left = T)
    
    # select only points outside polygons
    occ_out_nat <- occ_in_nat %>%
      filter(is.na(control)) %>%
      mutate(new_key = k, 
             new_species = unique(sp_info$new_species[sp_info$new_key==k])) %>%
      select(new_species, new_key, countryCode, geometry)
    
    # save in final data for next analysis
    # some species might have no record outside native range
    # thus save output only if at least one row in occ_out_nat
    
    if(nrow(occ_out_nat)>0){
      saveRDS(occ_out_nat, paste0(out_path, "RISK_15_alien_occ_IUCN_spk_", k))
    }
    print(k)
  }
  
}

length(list.files(out_path))


##### 3. Native country list from IUCN (n = 37 sp) ####

# file list native country lists
nat_path_iucn <- "Output/Native_exotic_range/Native_IUCN/"
nat_files_iucn <- list.files(nat_path_iucn)

sp_key_nat_countries <- sp_info %>% 
  filter(nat_exo_info %in% c("NAT_COUNTRY_IUCN")) %>%
  distinct(new_key) %>% pull(new_key)

# extract exotic occ for each sp in the list
for (k in sp_key_nat_countries){
  #k = sp_key_nat_countries[1]
  nat_k <- nat_files_iucn[grepl(paste0("spk_", k), nat_files_iucn)]
  out_k <- out_files[grepl(paste0("spk_", k), out_files)]
  print(k)
  
  if(is_empty(nat_k)){ next }
  if(is_empty(out_k)) {
    nat_list <- readRDS(paste0(nat_path_iucn, nat_k))
    
    occ_k <- occ_files[grepl(paste0("spk_", k), occ_files)]
    occ <- readRDS(paste0(occ_path, occ_k))
    
    # join dataframes by countryCode
    occ_in_nat <- left_join(occ, 
                            nat_list %>% rename(countryCode = code), 
                            by="countryCode")
    
    # select only points outside polygons
    occ_out_nat <- occ_in_nat %>%
      filter(is.na(binomial)) %>%
      mutate(new_key = k, 
             new_species = unique(sp_info$new_species[sp_info$new_key==k])) %>%
      # transform into sf format
      mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
      st_as_sf(
        coords = c("LONG", "LAT"),
        agr = "constant",
        crs = "+proj=longlat +datum=WGS84",
        stringsAsFactors = FALSE,
        remove = TRUE) %>%
      select(new_species, new_key, countryCode, geometry)
    
    # save in final data for next analysis
    # some species might have no record outside native range
    # thus save output only if at least one row in occ_out_nat
    
    if(nrow(occ_out_nat)>0){
      saveRDS(occ_out_nat, paste0(out_path, "RISK_15_alien_occ_IUCN_spk_", k))
    }
    print(k)
  }
  
}

length(list.files(out_path))


##### 4. Exotic country list from GRIIS (n = 99 sp) ####

# file list exotic country lists
exo_path <- "Output/Native_exotic_range/Exotic_GRIIS/"
exo_files <- list.files(exo_path)

sp_key_exo_countries <- sp_info %>% 
  filter(nat_exo_info %in% c("EXO_COUNTRY_GRIIS")) %>%
  distinct(new_key) %>% pull(new_key)


# extract exotic occ for each sp in the list
for (k in sp_key_exo_countries){
  #k = sp_key_exo_countries[50]
  exo_k <- exo_files[grepl(paste0("spk_", k), exo_files)]
  out_k <- out_files[grepl(paste0("spk_", k), out_files)]
  print(k)
  
  if(is_empty(exo_k)){ next }
  if(is_empty(out_k)) {
    exo_list <- readRDS(paste0(exo_path, exo_k))
    
    occ_k <- occ_files[grepl(paste0("spk_", k), occ_files)]
    occ <- readRDS(paste0(occ_path, occ_k))
    
    # join dataframes by countryCode
    occ_with_exo <- left_join(occ, exo_list, by="countryCode")
    
    # select only points outside polygons
    occ_in_exo <- occ_with_exo %>%
      filter(!is.na(new_species)) %>%
      mutate(new_key = k, 
             new_species = unique(sp_info$new_species[sp_info$new_key==k])) %>%
      # transform into sf format
      mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
      st_as_sf(
        coords = c("LONG", "LAT"),
        agr = "constant",
        crs = "+proj=longlat +datum=WGS84",
        stringsAsFactors = FALSE,
        remove = TRUE) %>%
      select(new_species, new_key, countryCode, geometry)

    
    # save in final data for next analysis
    # some species might have no record inside exotic range
    # thus save output only if at least one row in occ_in_exo
    
    if(nrow(occ_in_exo)>0){
      saveRDS(occ_in_exo, paste0(out_path, "RISK_15_alien_occ_IUCN_spk_", k))
    }
    print(k)
    print(which(sp_key_exo_countries==k))
  }
  
}

length(list.files(out_path))








######### BROUILLON ##########


nrow(test_occ)
test_occ <- readRDS("Data/True_exposure_alien_species/RISK_15_alien_occ_IUCN_spk_10830355")

test_nat <- readRDS("Output/Native_exotic_range/Native_IUCN/RISK_12_native_range_IUCN_spk_10830355")

ggplot() +
  geom_sf(data = test_nat) +
  geom_sf(data = test_occ)

test_occ <- readRDS("Data/True_exposure_alien_species/RISK_15_alien_occ_IUCN_spk_2439838")
test_nat <- readRDS("Output/Native_exotic_range/Native_IUCN/RISK_12_native_range_IUCN_spk_2439838")

ggplot() +
  geom_sf(data = world) +
  geom_sf(data = test_nat, fill = "firebrick") +
  geom_sf(data = test_occ)


test_occ <- readRDS("Data/True_exposure_alien_species/RISK_15_alien_occ_IUCN_spk_9752149")
test_nat <- readRDS("Output/Native_exotic_range/Native_IUCN/RISK_12_native_range_IUCN_spk_9752149")

library(spData)
data("world")


ggplot() +
  geom_sf(data = world)+
  geom_sf(data = test_nat, fill = "firebrick")+
  geom_sf(data = test_occ)





#######################
# brouillon

sp_with_pbm <- c("2489450", "5218911")
# try 
st_is_valid(st_make_valid(st_set_precision(nat_range, 1e5)))

which(sp_key=="2489450")
sp_key[77]
#####
# to remove geometry
occ %>%
  dplyr::mutate(LONG = sf::st_coordinates(.)[,1],
                LAT = sf::st_coordinates(.)[,2]) %>%
  st_set_geometry(NULL)