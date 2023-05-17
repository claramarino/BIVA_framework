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
world <- ne_countries(scale = "large", type = "countries", returnclass = "sf")


# attribute a country code (missing mostly for Namibia)
na_in_country <- st_intersection(st_make_valid(na_codes),
                                 st_make_valid(world)) 
na_in_country_ok <- na_in_country %>%
  mutate(countryCode = iso_a2) %>% 
  dplyr::select(new_species:coordinateUncertaintyInMeters) %>%
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
  left_join(sp_with_corresp, 
            sp_info %>% rename(ias_sp_name = original_sciname),
            relationship = "many-to-many"),
  left_join(sp_with_corresp, sp_info %>% rename(ias_sp_name = spe_lower),
            relationship = "many-to-many")) %>% 
  distinct(new_key, new_species, countryCode) %>%
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

# save outputs
saveRDS(comp_110, "Output/Exposure/RISK_22_completeness_GRIIS_DB_110km")
saveRDS(comp_55, "Output/Exposure/RISK_22_completeness_GRIIS_DB_55km")


##### Add sampling effort from GBIF #####

# utiliser la densité d'occurrences générée par meyer
# résolution = 100km

grid100km <- readRDS(paste0("Output/RISK_32_grid_110km"))

meyer = "Z:/THESE/5_Data/Distribution_spatiale/GBIF_data_bias_Meyer/"
list.files(meyer)

dl_meyer <- list(
  amph = read.csv2(paste0(meyer, "Meyer_etal_data_amph_110km.csv")) %>% 
    dplyr::select(Longitude, Latitude, RecordDensity) %>%
    rename(RecDens_amph = RecordDensity),
  mam = read.csv2(paste0(meyer, "Meyer_etal_data_mamm_110km.csv")) %>% 
    dplyr::select(Longitude, Latitude, RecordDensity) %>%
    rename(RecDens_mam = RecordDensity),
  bird = read.csv2(paste0(meyer, "Meyer_etal_data_aves_110km.csv")) %>% 
    dplyr::select(Longitude, Latitude, RecordDensity) %>%
    rename(RecDens_bird = RecordDensity)
) %>%
  reduce(full_join, by=c("Longitude", "Latitude")) %>%
  replace_na(list(RecDens_mam = 0, RecDens_amph = 0)) %>%
  dplyr::mutate(dens_all = RecDens_mam + RecDens_amph + RecDens_bird)

hist(log(dl_meyer$dens_all))

dens_sf <- sf::st_as_sf(dl_meyer,
                        coords = c("Longitude", "Latitude"),
                        crs = "+proj=longlat +datum=WGS84")

dens_sf_buff <- st_buffer(dens_sf, dist = 20000)

dens_sf_moll <- sf::st_transform(dens_sf_buff, 
                                 crs = raster::crs("+proj=moll")) %>%
  st_make_valid()

raster::crs(grid100km) == raster::crs(dens_sf_moll)

dens_grid <- st_intersection(grid100km, 
                             dens_sf_moll)

dens_cell <- dens_grid %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(dens = mean(dens_all)) %>%
  mutate(log_dens = log(dens+1))


dens_final <- left_join(grid100km, dens_cell)
hist(dens_final$dens)
hist(dens_final$log_dens)

ggplot(dens_final, aes(fill = log_dens))+
  geom_sf(color=NA)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar")



# Problem with some cells in the middle of continent with NA
# fill with mean value of neighbor cells
# how many NA ?
sum(is.na(dens_final$dens)) #3744

# for each NA, if cell touches other cells with value
# take the mean value of neighbors
# else, stay NA

# get cell neighbors
neighbors <- st_touches(dens_final)

# init: faster without geometry
dens_check <- dens_final %>% st_drop_geometry()

for (i in 1:nrow(dens_check)){
  if(is.na(dens_check$dens[i])){
    nbg <- neighbors[[i]]
    val <- dens_check$dens[nbg]
    if(sum(is.na(val)) < length(val)){
      dens_check$dens[i] <- min(na.omit(val))
    }
  }
}

sum(is.na(dens_final$dens))
sum(is.na(dens_check$dens)) # 1560


dens_final <- left_join(grid100km, dens_check) %>%
  mutate(log_dens = log(dens+1))

sum(is.na(dens_final$log_dens))

ggplot(dens_final, aes(fill = log_dens))+
  geom_sf(color=NA)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar")


dens_final <- replace_na(dens_final, list(log_dens = 0, dens = 0))

saveRDS(dens_final, "Output/Exposure/RISK_22_sampling_effort_gbif_meyer_110km")



############## Compute final completion Ce ###################

# Completion per country related to IAS (C)
# combined with sampling effort (e)

res = "110" # 110 or 55 km

se <- readRDS(paste0("Output/Exposure/RISK_22_sampling_effort_gbif_meyer_", res, "km"))
comp <- readRDS(paste0("Output/Exposure/RISK_22_completeness_GRIIS_DB_", res, "km"))

poly_iso <- st_read("Z:/THESE/5_Data/Distribution_spatiale/Countries_ISO3166_boundaries/iso3166-1-boundaries.shp")
str(poly_iso)

# associate grid with admin units
se_wgs84 <- st_transform(se, crs = st_crs(poly_iso)) %>%
  st_make_valid()
poly_iso_moll <- poly_iso %>%
  st_make_valid() %>%
  st_transform(crs = st_crs(se)) %>%
  st_make_valid()

ggplot(se_wgs84) + geom_sf()

se_iso <- st_intersection(se_wgs84, poly_iso)

se_iso <- left_join(se %>% dplyr::select(grid_id), 
                    se_iso %>% st_drop_geometry())

# final completeness per admin unit
tot <- comp %>% group_by(countryCode) %>%
  summarise(tot_sp_count = sum(nb_sp))

prop_country <- left_join(comp, tot) %>%
  mutate(prop = nb_sp/tot_sp_count) %>% 
  dplyr::select(countryCode, tot_sp_count, cate, prop) %>%
  pivot_wider(names_from = cate, values_from = prop) %>%
  filter(!is.na(tot_sp_count)) %>% 
  replace_na(list(`In our DB & GRIIS` = 0, `Only our DB` = 0)) %>%
  mutate(DB_GRIIS_or_not = `In our DB & GRIIS` + `Only our DB`) %>%
  mutate(None_in_GRIIS = if_else(`Only our DB` == 1, "None_in_GRIIS","Some_in_GRIIS")) %>%
  filter(!is.na(countryCode))


# associtae completeness and sampling effort

se_comp <- left_join(se_iso %>% rename(countryCode = ISO3166.1), prop_country)
sum(is.na(se_comp$log_dens))
colnames(se_comp)

se_comp_glob <- se_comp %>%
  st_drop_geometry() %>%
  group_by(grid_id, log_dens) %>%
  summarise(#log_dens=mean(log_dens),
            #DB_only = mean(`Only our DB`),
            #GRIIS_only = mean(`Only GRIIS`),
            DB_and_GRIIS = mean(`In our DB & GRIIS`),
            DB_GRIIS_or_not = mean(DB_GRIIS_or_not)) %>%
  filter(!is.na(DB_and_GRIIS)) %>%
  mutate(norm_dens = log_dens/max(se$log_dens)) %>%
  mutate(comp_se_sum = DB_and_GRIIS + norm_dens,
         comp_se_product = DB_and_GRIIS*norm_dens)

summary(se_comp_glob)


df_glob_plot <- left_join(se %>% dplyr::select(grid_id), se_comp_glob) %>%
  replace_na(list(comp_se_product = 0, comp_se_sum = 0))


# save final output

saveRDS(df_glob_plot,paste0("Output/Exposure/RISK_22_Comp_smpl_eff_", res, "km"))


#####################################

# brouillon => comparaison entre sum of components & product of components
# faire suppl ?


summary(df_glob_plot)

ggplot(df_glob_plot)+
  geom_sf(aes(fill = norm_dens), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)

ggplot(df_glob_plot)+
  geom_sf(aes(fill = DB_and_GRIIS), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)

ggplot(df_glob_plot)+
  geom_sf(aes(fill = comp_se_sum), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)

ggplot(df_glob_plot)+
  geom_sf(aes(fill = comp_se_product), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)

hist(df_glob_plot$comp_se_product)


# test cartogram 
library(cartogram)

df_5deg <- df_glob_plot %>% 
  st_make_grid(cellsize = 550000, square = F) %>% 
  st_sf()
ggplot(df_5deg)+geom_sf()

df_glob_5deg <- st_intersection(df_glob_plot, 
                                df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
  st_drop_geometry() %>%
  group_by(id_5deg) %>%
  summarise(comp_sum = mean(na.omit(comp_se_sum)),
            comp_prod = mean(na.omit(comp_se_product)))

df_glob_5deg_sf <- left_join(df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg)),
                             df_glob_5deg)
ggplot(df_glob_5deg_sf)+geom_sf(aes(fill=comp_prod))

st_crs(df_glob_5deg_sf)==st_crs(df_glob_plot)

cartog_prod <- cartogram_ncont(df_glob_plot,
                 weight = "comp_se_product")

# ne fonctionne pas en continu
cartog <- cartogram_cont(df_glob_5deg_sf[1:200,]%>% mutate(comp_prod=100000000*comp_prod),
                         weight = "comp_prod", itermax = 15)

ggplot(cartog_prod)+
  geom_sf()


ggplot(cartog_prod)+
  geom_sf(aes(fill=comp_se_product, color = comp_se_product))

hist(cartog_prod$comp_prod)
ggplot()+
  geom_sf(data = df_glob_5deg_sf[1:200,], color = "blue", alpha = .5)+
  geom_sf(data = cartog, color = "red")


# prop in griis & our db
ggplot(poly_prop) + geom_sf(aes(fill=`In our DB & GRIIS`), color = NA)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar") +
  geom_sf_text(aes(label=tot_sp_count), color ="black", size =3)

# prop in our db (GRIIS or not)
#install.packages("ggpattern")
library(ggpattern)

ggplot(poly_prop) + 
  geom_sf_pattern(aes(fill=DB_GRIIS_or_not, pattern = None_in_GRIIS), color = NA,
                  pattern_fill = "white", pattern_angle = 45, 
                  pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6)+
  scale_fill_gradient(high = "green", low = "red", guide = "colorbar") +
  scale_pattern_manual(values = c(None_in_GRIIS = "stripe", Some_in_GRIIS = "none"))+
  geom_sf_text(aes(label=tot_sp_count), color ="black", size =3)









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





#################### BROUILLON ################ 

# kernel smooth à faire autour des points ?
# générer un buffer autour de chaque point pour éviter les problème de non filling
# de certaines cellules ??

ggplot() +
  geom_sf(data = dens_sf_buff[1:3,])+
  geom_sf(data = dens_sf[1:3,])

ggplot() +
  geom_sf(data = dens_sf_buff, aes(fill = log(dens_all)))


dens_sf_buff <- st_buffer(dens_sf_moll, dist = 20000)
dens_sf_buff

zone_xy <- grid100km %>% 
  st_make_valid() %>%
  st_transform(crs = "+proj=longlat +datum=WGS84") %>%
  st_centroid() %>%
  st_coordinates() %>% 
  as_tibble() %>% 
  dplyr::select(x = X, y = Y)


to_smooth <- dens_sf %>% 
  cbind(., st_coordinates(.)) %>%
  st_set_geometry(NULL) %>% 
  dplyr::select(x = X, y = Y, dens_all) 

ks <- btb::btb_smooth(pts = to_smooth,
                      sEPSG = raster::crs(dens_sf),
                      iCellSize = 1L,
                      iBandwidth = 20L,
                      vQuantiles = NULL) %>%
  st_set_crs(raster::crs("+proj=longlat +datum=WGS84"))%>%
  st_transform(raster::crs("+proj=moll"))

ggplot(ks)+ geom_sf(aes(fill = dens_all), color=NA)



install.packages("btb")
library(btb)



