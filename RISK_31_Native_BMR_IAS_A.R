# Retrieve native range for all species in iucn that are IAS-A
# + filter species from ou ias_x_native list


rm(list=ls())

library(tidyverse)
library(stringr)
library(readr)
library(sp)
library(sf)
library(rredlist)
library(raster)
library(rnaturalearth) # devtools::install_github("ropensci/rnaturalearthhires")
library(pbmcapply)

#### Load data and extract ias A/T info ####

# load threats
native_iasa_x_ias <- readRDS("characterize_ias_yan/Data/00_all_native_ias_a")
length(unique(native_iasa_x_ias$scientificName))
table(native_iasa_x_ias %>% distinct(scientificName, className) %>% pull(className))

# compute T, NT, EX-EW col
native_iasa <- native_iasa_x_ias %>%
  distinct(scientificName, className, redlistCategory) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Endangered", "Critically Endangered", "Vulnerable"),
    "T", "NT_DD")) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Extinct", "Extinct in the Wild"), 
    "EX_EW", threatened))


# add a column for species in ias_x_native with specified ias
native_spe_ias <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok") %>%
  distinct(scientificName) %>% pull(scientificName)

native_iasa <- native_iasa %>%
  mutate(ias_modelled = if_else(scientificName %in% native_spe_ias, "YES","NO"))




##### Polygons of IAS A #####

in_fold = "Z:/THESE/5_Data/Distribution_spatiale/5_IUCN_MAJ_2022_07_AMPH_REPT_MAM/"
out_fold = "Output/Sensitivity/Polygons_native_IAS_A"

#### MAMMALS

# load native ranges from IUCN
mam_all <- st_read(paste0(in_fold, "MAMMALS.shp"))

mam_nat <- mam_all %>% 
  dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% native_iasa$scientificName)%>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry) %>%
  dplyr::filter(terrestial == "true" | freshwater == "true")

head(mam_nat)

length(unique(mam_nat$binomial))
table(native_iasa$className)

native_iasa %>% 
  dplyr::filter(className =="MAMMALIA") %>%
  dplyr::filter(!(scientificName %in% unique(mam_nat$binomial)))
# all missing species are EX_EW

# save native polygons all ias_a mammals
mam_nat_iasa <- left_join(mam_nat, native_iasa %>%
                            select(scientificName, threatened, ias_modelled) %>%
                            rename(binomial = scientificName),
                          by = "binomial")
saveRDS(mam_nat_iasa, paste0(out_fold, "/RISK_31_Polygon_nat_ias_a_MAM"))


#### REPTILES

# load native ranges from IUCN
rept_all <- st_read(paste0(in_fold, "REPTILES.shp"))

rept_nat <- rept_all %>% 
  dplyr::filter(origin %in% c(1, 2, 5, 6)) %>%
  dplyr::filter(binomial %in% native_iasa$scientificName)%>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry)%>%
  dplyr::filter(terrestial == "true" | freshwater == "true")

length(unique(rept_nat$binomial)) # 760
table(native_iasa$className) # 834

miss_sp <- native_iasa %>% 
  dplyr::filter(className =="REPTILIA") %>%
  dplyr::filter(!(scientificName %in% unique(rept_nat$binomial))) %>% 
  pull(scientificName)
# some missing sp => even not in IUCN MAJ, check in GARD?
rm(rept_all)
rept_gard <- st_read("Z:/THESE/5_Data/Distribution_spatiale/REPTILE_Distrib_GARD.1_dissolved_ranges/modeled_reptiles.shp")
head(rept_gard)

gard_filt <- rept_gard %>% 
  dplyr::filter(Binomial %in% miss_sp) %>%
  dplyr::rename(binomial = Binomial) %>%
  dplyr::select(binomial, geometry)
length(unique(gard_filt$binomial))

rept_nat_all <- bind_rows(rept_nat, gard_filt)

# save native polygons all ias_a reptiles
rept_nat_iasa <- left_join(rept_nat_all, native_iasa %>%
                            select(scientificName, threatened, ias_modelled) %>%
                            rename(binomial = scientificName),
                          by = "binomial")

saveRDS(rept_nat_iasa, paste0(out_fold, "/RISK_31_Polygon_nat_ias_a_REPT"))


#### BIRDS 

# load native ranges from Birdlife
path_birds <- "Z:/THESE/5_Data/Distribution_spatiale/1_MAJ_BIRDS_Birdlife/"

# need to correct names pbm
good_names <- c("Antigone antigone" = "Grus antigone",
                "Geomalia heinrichi" = "Zoothera heinrichi",
                "Chatarrhaea longirostris" = "Argya longirostris",
                "Chaetornis striata" = "Schoenicola striatus")

b_names_ok <- c(native_iasa$scientificName, names(good_names))

b1 <- st_read(paste0(path_birds, "Birds_A_2600_repare.shp"))
head(b1)
b1_filt <- b1 %>% dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% b_names_ok)%>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b1)
b2 <- st_read(paste0(path_birds, "Birds_A_3400_repare.shp"))
b2_filt <- b2 %>% dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% b_names_ok)%>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b2)
b3 <- st_read(paste0(path_birds, "Birds_B_5800_repare.shp"))
b3_filt <- b3 %>% dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% b_names_ok)%>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b3)
b4 <- st_read(paste0(path_birds, "Birds_C_5766_repare.shp"))
b4_filt <- b4 %>% dplyr::filter(origin %in% c(1,2)) %>%
  dplyr::filter(binomial %in% b_names_ok)%>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b4)

# bind the four files together
bird_nat <- bind_rows(b1_filt, b2_filt, b3_filt, b4_filt)

bird_nat$binomial[bird_nat$binomial %in% names(good_names)] <-
  good_names[bird_nat$binomial[bird_nat$binomial %in% names(good_names)]]


length(unique(bird_nat$binomial)) # 876
table(native_iasa$className) # 877

native_iasa %>% 
  dplyr::filter(className =="AVES") %>%
  dplyr::filter(!(scientificName %in% unique(bird_nat$binomial)))

# save native polygons all ias_a birds
bird_nat_iasa <- left_join(bird_nat, native_iasa %>%
                             select(scientificName, threatened, ias_modelled) %>%
                             rename(binomial = scientificName),
                           by = "binomial")

saveRDS(bird_nat_iasa, paste0(out_fold, "/RISK_31_Polygon_nat_ias_a_BIRD"))




#### Extract overlapping cells from raster ####

# path for loading polygons
path_poly = "~/predict_vulnerability/Polygons_native_IAS_A/"

# define raster
# high resolution for keeping most islands
# initialize raster grids
# resolution 0.1
raster0.1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    resolution = 0.1, vals=NULL)
raster0.1[] <- 0
# resolution 1
raster1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                  resolution = 1, vals=NULL)
raster1[] <- 0


# world map to filter out oceanic ranges
worldMap <- ne_countries(scale = "large", type = "countries", returnclass = "sf")
# filter directly rasters
raster0.1 <- mask(raster0.1, worldMap)
raster1 <- mask(raster1, worldMap)
  

# define extract_cells function ----------

# INPUT : one polygon from multiple poly object + one raster grid with values = 1 
# OUTPUT : a df with x,y coordinates of all grid cells that are included in the poly
# and that are included in a land range (exclude oceanic ranges)
extract_cells <- function(poly, rast){
  # from a raster containing 1 for each cell
  # keep cells that are masked by the polygon = cells from native range 
  masked_rast <- mask(rast, poly)
  
  # for saving memory, save a df containing all cells with values
  df_val <- as.data.frame(masked_rast, xy = TRUE) %>%
    # remove empty cells
    dplyr::filter(!is.na(layer)) %>%
    dplyr::mutate(binomial = poly$binomial)
  
  return(df_val)
}

# folder for saving outputs
out_path <- "~/predict_vulnerability/"

# Computation ----------------

# apply the function to each line of Large spatial polygon df object
# separate per class


#### mammals ####

mam_ias_as <- readRDS(paste0(path_poly, "RISK_31_Polygon_nat_ias_a_MAM"))
mam_ias_as_ok <- mam_ias_as %>%
  dplyr::select(binomial, presence:seasonal, threatened:geometry) %>%
  dplyr::mutate(objectid = 1:nrow(mam_ias_as))
rm(mam_ias_as)

mam_IDs <- as.list(mam_ias_as_ok$objectid)

# for raster at 0.1 degres
df_mam_0.1_mclpply <- pbmclapply(
  mam_IDs,
  function(x){
    poly <- mam_ias_as_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 15)
saveRDS(df_mam_0.1_mclpply, paste0(out_path, "RISK_31_cells_nat_ias_a_MAM_r01"))

# for raster at 1 degres 
df_mam_1_mclpply <- pbmclapply(
  mam_IDs,
  function(x){
    poly <- mam_ias_as_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 15)

saveRDS(df_mam_1_mclpply, paste0(out_path, "RISK_31_cells_nat_ias_a_MAM_r1"))


#### reptiles ####

rept_ias_as <- readRDS(paste0(path_poly, "RISK_31_Polygon_nat_ias_a_REPT"))
rept_ias_as_ok <- rept_ias_as %>%
  dplyr::select(binomial, presence:seasonal, threatened:geometry) %>%
  dplyr::mutate(objectid = 1:nrow(rept_ias_as))
rm(rept_ias_as)

rept_IDs <- as.list(rept_ias_as_ok$objectid)

df_rept_0.1_mclpply <- pbmclapply(
  rept_IDs,
  function(x){
    poly <- subset(rept_ias_as_ok, objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 15)
saveRDS(df_rept_0.1_mclpply, paste0(out_path, "RISK_31_cells_nat_ias_a_REPT_r01"))

df_rept_1_mclpply <- pbmclapply(
  rept_IDs,
  function(x){
    poly <- subset(rept_ias_as_ok, objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 15)
saveRDS(df_rept_1_mclpply, paste0(out_path, "RISK_31_cells_nat_ias_a_REPT_r1"))

#### birds ####

birds_ias_as <- readRDS(paste0(path_poly, "RISK_31_Polygon_nat_ias_a_BIRD"))
bird_ias_as_ok <- birds_ias_as %>%
  dplyr::select(binomial, presence:seasonal, threatened:geometry) %>%
  dplyr::mutate(objectid = 1:nrow(birds_ias_as))
rm(birds_ias_as)

bird_IDs <- as.list(bird_ias_as_ok$objectid)

df_bird_0.1_pblpply <- pbmclapply(
  bird_IDs,
  function(x){
    poly <- subset(bird_ias_as_ok, objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 10)
saveRDS(df_bird_0.1_pblpply, paste0(out_path, "RISK_31_cells_nat_ias_a_BIRD_r01"))
# Elapsed 27:17

df_bird_1_pblpply <- pbmclapply(
  bird_IDs,
  function(x){
    poly <- subset(bird_ias_as_ok, objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 5)
saveRDS(df_bird_1_pblpply, paste0(out_path, "RISK_31_cells_nat_ias_a_BIRD_r1"))



#############################

# Message d'avis :
# Dans mclapply(rept_IDs, function(x) { :
#   le coeur programmé 10 n’a pas renvoyé de résultat, 
# toutes les valeurs pour cette tâche seront affectées
# 
# for (i in 1:length(df_rept_01_mclpply)){
#   if(is.null(df_rept_01_mclpply[[i]])){
#     poly <- subset(rept_ias_as, OBJECTID == i)
#     df_rept_01_mclpply[[i]] <- extract_cells(poly, raster0.1)
#     print(i)
#   }
# }

mam_df <- bind_rows(df_bird_0.1_pblpply) %>% distinct()

m_all_agg <- mam_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())


SR_m_as<- ggplot(data = m_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_m_as

