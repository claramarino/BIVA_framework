# find exotic/native range of species not in GRIIS and not in IUCN

rm(list=ls())

library(tidyverse)
library(sf)
library(sp)
library(readr)
library(rnaturalearth)


sp_info <- readRDS("Output/Native_exotic_range/RISK_13_ias_list_with_occ_IUCN_GRIIS")
sp_info <- sp_info %>% mutate(nat_exo_info = if_else(
  nat_exo_info=="NAT_COUNTRY_LIST", "NAT_COUNTRY_IUCN", nat_exo_info))

missing_sp <- sort(sp_info %>% filter(is.na(nat_exo_info)) %>% 
                     pull(new_species))

missing_sp


##### find range in CABI ##########

# "Neustanthus phaseoloides" = "Pueraria montana var. lobata" # also known as kudzu
# Miconia hirta = Clidemia hirta 

sp_in_cabi <- c("Batis maritima" = "Batis maritima", 
                "Miconia hirta" = "Clidemia hirta",
                "Neustanthus phaseoloides" = "Pueraria montana var. lobata",
                "Schistocerca nitens" = "Schistocerca nitens", 
                "Vespula pensylvanica" = "Vespula pensylvanica")

cabi_fold <- "Z:/THESE/5_Data/Distribution_spatiale/Alien_species/CABI_ISC_Manual_download/"
list.files(cabi_fold)


# get countries or state polygons from natural earth data

# polygons form countries
poly_countr <- rnaturalearth::ne_countries(returnclass = "sf", type = "map_units") %>%
  select(name, geometry)
poly_countr <- st_make_valid(poly_countr)

# polygons from states
poly_state <- rnaturalearth::ne_states(returnclass = "sf") %>%
  select(name_en, geometry)
poly_state <- st_make_valid(poly_state)


for(sp in sp_in_cabi){
  
  #sp=sp_in_cabi[2]
  
  cabi_range <- read_csv(paste0(cabi_fold, sp, ".csv"), skip = 1)
  native_cabi <- cabi_range %>% filter(Origin=="Native") %>%
    mutate(species = sp)
  
  
  # select all real countries (no tiret)
  countr_cab <- native_cabi[!grepl("- ", native_cabi$Region),]
  # select all states (with tiret) and remove tiret for good writing
  state_cab <- native_cabi[grepl("- ", native_cabi$Region),]
  
  
  # for real countries, export country poly with points in it
  # convert occ data frame into sf object
  occ_sf <- countr_cab %>%
    mutate_at(vars(Longitude, Latitude), as.numeric) %>%   # coordinates must be numeric
    st_as_sf(
      coords = c("Longitude", "Latitude"),
      agr = "constant",
      crs = "+proj=longlat +datum=WGS84",
      stringsAsFactors = FALSE,
      remove = TRUE)
  # join points and polygons 
  occ_in_count <- st_join(occ_sf, poly_countr, join = st_within, left = T)
  # search for "countries" without polygon
  no_count <- occ_in_count %>% filter(is.na(name)) %>% select(-name)
  
  
  # for provinces, export state poly with points in it
  # convert occ data frame into sf object
  occ_sf_state <- state_cab %>%
    mutate_at(vars(Longitude, Latitude), as.numeric) %>%   # coordinates must be numeric
    st_as_sf(
      coords = c("Longitude", "Latitude"),
      agr = "constant",
      crs = "+proj=longlat +datum=WGS84",
      stringsAsFactors = FALSE,
      remove = TRUE)
  
  occ_sf_state_no_count <- bind_rows(occ_sf_state, no_count)
  # join points and polygons 
  occ_in_state <- st_join(occ_sf_state_no_count, poly_state, join = st_within, left = T)
  
  final_range <- bind_rows(
    poly_countr %>% filter(name %in% occ_in_count$name),
    poly_state %>% filter(name_en %in% occ_in_state$name_en) %>%
      rename(name = name_en)
  )
  
  spk = unique(sp_info$new_key[sp_info$new_species==names(sp_in_cabi[which(sp_in_cabi==sp)])])
  
  saveRDS(final_range, 
          paste0("Output/Native_exotic_range/Other_database/",
                 "/RISK_14_native_state_CABI_spk_", spk))
  
}

# save info regarding data source
sp_info$nat_exo_info[sp_info$new_species %in% names(sp_in_cabi)] <- "NAT_COUNTRY_CABI"


# Cas par cas pour les sp hors CABI

################
# Anolis wattsii = Anolis wattsi in iucn ?
rept_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/2_MAJ_REPTILES_IUCN_2022/REPTILES.shp")
# select only the 3 species and native range 
poly_anolis <- subset(rept_all, binomial == "Anolis wattsi" &
                           origin %in% c(1,2))
key_anolis <- sp_info$new_key[sp_info$new_species=="Anolis wattsii"]
saveRDS(poly_anolis, paste0("Output/Native_exotic_range/Native_IUCN",
                        "/RISK_14_native_range_IUCN_spk_", key_anolis))
# save info regarding data source
sp_info$nat_exo_info[sp_info$new_key==key_anolis] <- "NAT_RANGE_POLY"

#################
# morelia imbricata = subspecies of morelia spilota
# use sp range from iucn for morelia spilota as native range
poly_morelia <- subset(rept_all, binomial == "Morelia spilota" &
                      origin %in% c(1,2))
key_morelia <- sp_info$new_key[sp_info$new_species=="Morelia imbricata"]
saveRDS(poly_morelia, paste0("Output/Native_exotic_range/Native_IUCN",
                        "/RISK_14_native_range_IUCN_spk_", key_morelia))
#save info of datasource
sp_info$nat_exo_info[sp_info$new_key==key_morelia] <- "NAT_RANGE_POLY"
rm(rept_all)

#################
# bos frontalis = bos gaurus in iucn, use polygon range map
mam_all <- st_read("Z:/THESE/5_Data/Distribution_spatiale/0_MAJ_IUCN_2020/Mam_repare.shp")
# select only the 3 species and native range 
poly_bos <- subset(mam_all, binomial == "Bos gaurus" &
                        origin %in% c(1,2))
key_bos <- sp_info$new_key[sp_info$new_species=="Bos frontalis"]
saveRDS(poly_bos, paste0("Output/Native_exotic_range/Native_IUCN",
                            "/RISK_14_native_range_IUCN_spk_", key_bos))
# save info regarding data source
sp_info$nat_exo_info[sp_info$new_key==key_bos] <- "NAT_RANGE_POLY"



#################
# lea floridensis
# native from florida # https://orthsoc.org/sina/131a.htm
# put florida as native country?

usa <- rnaturalearth::ne_states(country =  "United States of America", returnclass = "sf")
florida <- usa %>% filter(name=='Florida')%>% 
  select(name_en, geometry) %>% rename(name = name_en)
plot(florida)

key_lea <- sp_info$new_key[sp_info$new_species=="Lea floridensis"]
saveRDS(florida, paste0("Output/Native_exotic_range/Other_database",
                            "/RISK_14_native_state_florida_spk_", key_lea))

sp_info$nat_exo_info[sp_info$new_key==key_lea] <- "NAT_COUNTRY_HAND"


#################
# ceodes umbellata = Pisonia umbellifera
# it is native to the Andaman Islands, Indonesia, Malaysia, 
# the Philippines, Thailand, Vietnam, China, Taiwan, Hawaii, 
# Africa and Madagascar and the states of New South Wales and Queensland in Australia.

australia <- rnaturalearth::ne_states(country = "australia", returnclass = "sf")
india <- rnaturalearth::ne_states(country = "india", returnclass = "sf")

ceodes_usa <- usa %>% filter(name=='Hawaii')
ceodes_aust <- australia %>% filter(name %in% c("New South Wales" , "Queensland"))
ceodes_india <- india %>% filter(name_en=="Andaman and Nicobar Islands")

ceodes_count <- rnaturalearth::ne_countries(
  country = c("Indonesia", "Malaysia", "Philippines", "Thailand", "Vietnam", "China",
              "Madagascar", "Taiwan", "Japan"), returnclass = "sf")

ceodes_tot <- bind_rows(
  bind_rows(ceodes_usa, ceodes_aust, ceodes_india) %>% 
    select(name_en, geometry) %>% rename(name = name_en),
  ceodes_count %>% select(name, geometry))

plot(ceodes_tot)

key_ceodes <- sp_info$new_key[sp_info$new_species=="Ceodes umbellata"]
saveRDS(ceodes_tot, paste0("Output/Native_exotic_range/Other_database",
                        "/RISK_14_native_countries_hand_spk_", key_ceodes))

sp_info$nat_exo_info[sp_info$new_key==key_ceodes] <- "NAT_COUNTRY_HAND"



############ SAVE INFO TABLE #############
table(sp_info %>% distinct(new_key, nat_exo_info) %>% pull(nat_exo_info))

saveRDS(sp_info, "Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")

