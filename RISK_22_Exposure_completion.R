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
# ggplot(na_codes) +
#   geom_sf(data = world)+
#   geom_sf()
# 
# st_crs(na_codes)
# st_crs(world)

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

# keep only the species in the 310 that have exotic occurrences
out_fold = "Output/Exposure/Exposure_raw/"
df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km_nofilter"))
sp_310 <- unique(df_all$new_species)

length(unique(sp_with_corresp$ias_sp_name))

sp_info_in_griis <- bind_rows(
  left_join(sp_with_corresp, sp_info %>% rename(ias_sp_name = original_sciname)),
  left_join(sp_with_corresp, sp_info %>% rename(ias_sp_name = spe_lower))
) %>% distinct(new_key, new_species, countryCode) %>%
  filter(new_species %in% sp_310) %>%
  mutate(InGriis="yes") %>%
  filter(!is.na(new_key)) %>% 
  # some country codes are wrong in GRIIS
  # keep only the two first letters
  mutate(countryCode = strtrim(countryCode, 2)) %>% distinct()

length(unique(sp_info_in_griis$new_key))

# 273 species for which we have a country information

str(sp_country)
length(unique(sp_country$new_species))
length(unique(sp_info_in_griis$new_species))

setdiff(unique(sp_info_in_griis$new_species), unique(sp_country$new_species))
setdiff(unique(sp_country$new_species), unique(sp_info_in_griis$new_species))


comp_nolimit <- full_join(sp_country, sp_info_in_griis) %>% distinct() %>%
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



comp_55 <- full_join(sp_country_55, sp_info_in_griis) %>% distinct() %>%
  # keep only the 284 sp for which we have info in GRIIS
  filter(new_species %in% sp_info_in_griis$new_species) %>%
  replace_na(list(InOurDB = "no", InGriis = "no")) %>%
  mutate(cate = case_when(
    InOurDB == "yes" & InGriis =="yes" ~ "In our DB & GRIIS",
    InOurDB == "no" & InGriis =="yes" ~ "Only GRIIS",
    InOurDB == "yes" & InGriis =="no" ~ "Only our DB"
  )) %>% group_by(countryCode, cate) %>%
  summarise(nb_sp = n())

tot55 <- comp_55 %>% group_by(countryCode) %>%
  summarise(tot_sp_count = sum(nb_sp))

prop_country55 <- left_join(comp_55, tot55)%>%
  mutate(prop = nb_sp/tot_sp_count)

comp_110 <- full_join(sp_country_110, sp_info_in_griis) %>% distinct() %>%
  # keep only the 284 sp for which we have info in GRIIS
  filter(new_species %in% sp_info_in_griis$new_species) %>%
  replace_na(list(InOurDB = "no", InGriis = "no")) %>%
  mutate(cate = case_when(
    InOurDB == "yes" & InGriis =="yes" ~ "In our DB & GRIIS",
    InOurDB == "no" & InGriis =="yes" ~ "Only GRIIS",
    InOurDB == "yes" & InGriis =="no" ~ "Only our DB"
  )) %>% group_by(countryCode, cate) %>%
  summarise(nb_sp = n())

tot110 <- comp_110 %>% group_by(countryCode) %>%
  summarise(tot_sp_count = sum(nb_sp))

prop_country110 <- left_join(comp_110, tot110)%>%
  mutate(prop = nb_sp/tot_sp_count)





############## Map countries completion ###################

poly_iso <- st_read("Z:/THESE/5_Data/Distribution_spatiale/Countries_ISO3166_boundaries/iso3166-1-boundaries.shp")
str(poly_iso)

#ggplot(poly_iso) + geom_sf()

poly_prop <- left_join(poly_iso %>% rename(countryCode = ISO3166.1), 
                      prop_country %>% select(countryCode, tot_sp_count, cate, prop) %>%
                        pivot_wider(names_from = cate, values_from = prop)) %>%
  filter(!is.na(tot_sp_count)) %>% 
  replace_na(list(`In our DB & GRIIS` = 0, `Only our DB` = 0)) %>%
  mutate(DB_GRIIS_or_not = `In our DB & GRIIS` + `Only our DB`) %>%
  mutate(None_in_GRIIS = if_else(`Only our DB` == 1, "None_in_GRIIS","Some_in_GRIIS"))

colnames(poly_prop)


poly_prop55 <- left_join(poly_iso %>% rename(countryCode = ISO3166.1), 
                         prop_country55 %>% select(countryCode, tot_sp_count, cate, prop) %>%
                           pivot_wider(names_from = cate, values_from = prop)) %>%
  filter(!is.na(tot_sp_count)) %>% 
  replace_na(list(`In our DB & GRIIS` = 0, `Only our DB` = 0)) %>%
  mutate(DB_GRIIS_or_not = `In our DB & GRIIS` + `Only our DB`) %>%
  mutate(None_in_GRIIS = if_else(`Only our DB` == 1, "None_in_GRIIS","Some_in_GRIIS"))


poly_prop110 <- left_join(poly_iso %>% rename(countryCode = ISO3166.1), 
                          prop_country110 %>% select(countryCode, tot_sp_count, cate, prop) %>%
                            pivot_wider(names_from = cate, values_from = prop)) %>%
  filter(!is.na(tot_sp_count)) %>% 
  replace_na(list(`In our DB & GRIIS` = 0, `Only our DB` = 0)) %>%
  mutate(DB_GRIIS_or_not = `In our DB & GRIIS` + `Only our DB`) %>%
  mutate(None_in_GRIIS = if_else(`Only our DB` == 1, "None_in_GRIIS","Some_in_GRIIS"))


# prop in griis & our db
ggplot(poly_prop) + geom_sf(aes(fill=`In our DB & GRIIS`), color = NA)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar") +
  geom_sf_text(aes(label=tot_sp_count), color ="black", size =3)

# prop in our db (GRIIS or not)
#install.packages("ggpattern")
library(ggpattern)

ggplot(poly_prop55) + 
  geom_sf_pattern(aes(fill=DB_GRIIS_or_not, pattern = None_in_GRIIS), color = NA,
                  pattern_fill = "white", pattern_angle = 45, 
                  pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar") +
  scale_pattern_manual(values = c(None_in_GRIIS = "stripe", Some_in_GRIIS = "none"))+
  geom_sf_text(aes(label=tot_sp_count), color ="black", size =3)


########################
# output folder
out_fold = "Output/Exposure/Exposure_raw/"

# df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km_nofilter"))
df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km"))
df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_55km"))



# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())

df_all_range <- left_join(df_all %>% dplyr::select(-new_key), 
                          tot_range,
                          by = "new_species")

# calculate cell related metrics
ias_all_agg <- df_all_range %>%
  group_by(grid_id) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_occ), # nb occ per cell (all ias sp)
    nb_occ_med = median(n_occ), # median nb occ per cell
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range)
  ) 

# Using Berhmann's cylindrical equal-area projection (CEA)
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load worldmap for spatial extent 
world <- rnaturalearth::countries110
# Convert into sf objects with CEA projection
study_area_sf <- sf::st_transform(sf::st_as_sf(world), crs=cea)
st_bbox(world)
st_bbox(study_area_sf)


# create grid cells 110 km
grid_110km = st_make_grid(
  x= study_area_sf, cellsize = 110000, what = "polygons", square = F) %>%
  st_sf() 
grid_110km = grid_110km %>%
  mutate(grid_id = 1:nrow(grid_110km))


ias_agg_grid110 <- st_as_sf(left_join(ias_all_agg, grid_110km))

# plot species richness
SR_ias<- ggplot(data = ias_agg_grid110) +
  geom_sf(aes(fill = SR_tot), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1) +
  geom_sf(data = study_area_sf, alpha = 0.1, fill = NA) +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$new_key))),
       x = "Longitude", y = "Latitude")
SR_ias

poly_prop110_cea <- sf::st_transform(poly_prop110, crs=cea) 

st_bbox(poly_prop110_cea)
st_bbox(grid_110km)

inter_grid_compl <- st_intersection(grid_110km, poly_prop110_cea)

str(inter_grid_compl)
length(unique(inter_grid_compl$grid_id))

inter_value <- inter_grid_compl %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(
    min_only_our_db = min(Only.our.DB),
    min_only_griis = min(Only.GRIIS),
    min_db_and_griis = min(In.our.DB...GRIIS),
    min_db_griisornot= min(DB_GRIIS_or_not),
    mean_only_our_db = mean(Only.our.DB),
    mean_only_griis = mean(Only.GRIIS),
    mean_db_and_griis = mean(In.our.DB...GRIIS),
    mean_db_griisornot= mean(DB_GRIIS_or_not)
  )

grid_val <- left_join(grid_110km, inter_value) %>% filter(!is.na(min_only_our_db))

ggplot(grid_val)+ geom_sf(aes(fill=min_db_griisornot), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)


hist(grid_val$min_db_and_griis)
ggplot(grid_val)+
  geom_histogram(aes(x = min_db_griisornot), fill = "red")+
  geom_histogram(aes(x = min_db_and_griis), alpha = .5)



######################################################################

# How ok our data are for climate SDM? probably OK tiede


# get all GBIF points related to each species (native and alien)
# count the number of cells they fall in for cliamte SDM

# filtering for different precisions?

rm(list=ls())

library(tidyverse)
library(raster)
library(stringr)
library(readr)
library(sf)

# species key
out_fold = "Output/Exposure/Exposure_raw/"
df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km_nofilter"))
sp_310 <- unique(df_all$new_key)

# file list occurrence points
occ_path <- "Output/Occurrences_clean_taxo_ok/"
occ_files <- list.files(occ_path)


#open chelsa raster

chelsa <- raster("Data/CHELSA_bio10_01.tif")
plot(chelsa)
str(chelsa)

chelsa_10arcm<- raster::aggregate(chelsa, 
                                 fact = 20)

chelsa_10arcm
plot(chelsa_10arcm)

df <- data.frame()

for (k in sp_310[50:310]){
  #k = "2498305"
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
  
  c_pts <- as.data.frame(extract(chelsa_10arcm, occ_sf, cellnumbers=TRUE))
  
  occ110 <- occ_sf %>% filter(coordinateUncertaintyInMeters < 55000)
  c110 <- as.data.frame(extract(chelsa_10arcm, occ110, cellnumbers=TRUE))
  
  occ55 <- occ_sf %>% filter(coordinateUncertaintyInMeters < 22500)
  c55 <- as.data.frame(extract(chelsa_10arcm, occ55, cellnumbers=TRUE))
  
  dfk = data.frame(sp_key = k,
                   n_pts_all = nrow(occ_sf),
                   n_pts_110 = nrow(occ110),
                   n_pts_55 = nrow(occ55),
                   n_cells_all = length(unique(c_pts$cells)),
                   n_cells_110 = length(unique(c110$cells)),
                   n_cells_55 = length(unique(c55$cells)))
  
  df <- bind_rows(df, dfk)
  print(k)
}
# bug for k=49

ggplot(df)+geom_histogram(aes(log(n_cells_110)))+
  geom_vline(xintercept = log(30))
nrow(df%>% filter(n_cells_110>30))

ggplot(df)+geom_histogram(aes(log(n_cells_55)))+
  geom_vline(xintercept = log(30))
nrow(df%>% filter(n_cells_55>30))


saveRDS(df, "Output/Exposure/RISK_22_completion_n_cells_climate")
