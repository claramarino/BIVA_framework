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
sp_info <- readRDS("Output/Native_exotic_range/RISK_11_ias_list_with_occ_IUCN_check")

sp_key <- unique(sp_info$new_key)


# file list occurrence points
occ_path <- "Output/Occurrences_clean_taxo_ok/"
occ_files <- list.files(occ_path)
# file list native ranges
nat_path <- "Output/Native_exotic_range/Native_IUCN/"
nat_files <- list.files(nat_path)

# file list output
# do only sp for which you haven't done it yet
out_path <- "Data/True_exposure_alien_species/"
# don't forget to clean folder if doing again the procedure
out_files <- list.files(out_path)
# out_files <- character()

for (k in sp_key){
  # k = "2498305"
  nat_k <- nat_files[grepl(paste0("spk_", k), nat_files)]
  out_k <- out_files[grepl(paste0("spk_", k), out_files)]
  print(k)
  
  if(is_empty(nat_k)){ next }
  if(is_empty(out_k)) {
    nat_range <- readRDS(paste0(nat_path, nat_k))
    
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
    
    # create a control column for removing occurrences outside poly
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




# open one file


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