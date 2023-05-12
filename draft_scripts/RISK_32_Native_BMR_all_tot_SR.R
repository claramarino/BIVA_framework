# Retrieve native range for all species in iucn 
# total species richness

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

##### Polygons of Native, not ex #####

in_fold = "Z:/THESE/5_Data/Distribution_spatiale/5_IUCN_MAJ_2022_07_AMPH_REPT_MAM/"
out_fold = "Output/Sensitivity/Polygons_native_all_bmr"

#### MAMMALS

# load native ranges from IUCN
mam_all <- st_read(paste0(in_fold, "MAMMALS.shp"))

mam_nat <- mam_all %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry) %>%
  dplyr::filter(terrestial == "true" | freshwater == "true") %>%
  dplyr::filter(category %in% c("CR","EN","VU","LC","DD","NT"))


saveRDS(mam_nat, paste0(out_fold, "/RISK_32_Polygon_nat_all_MAM"))


#### REPTILES

# load native ranges from IUCN
rept_all <- st_read(paste0(in_fold, "REPTILES.shp"))

rept_nat <- rept_all %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry) %>%
  dplyr::filter(terrestial == "true" | freshwater == "true") %>%
  dplyr::filter(category %in% c("CR","EN","VU","LC","DD","NT"))

saveRDS(rept_nat, paste0(out_fold, "/RISK_32_Polygon_nat_all_REPT"))
head(rept_nat)

#### BIRDS 

# load native ranges from Birdlife
path_birds <- "Z:/THESE/5_Data/Distribution_spatiale/1_MAJ_BIRDS_Birdlife/"

b1 <- st_read(paste0(path_birds, "Birds_A_2600_repare.shp"))
b1_filt <- b1 %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b1)
b2 <- st_read(paste0(path_birds, "Birds_A_3400_repare.shp"))
b2_filt <- b2  %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b2)
b3 <- st_read(paste0(path_birds, "Birds_B_5800_repare.shp"))
b3_filt <- b3  %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b3)
b4 <- st_read(paste0(path_birds, "Birds_C_5766_repare.shp"))
b4_filt <- b4 %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(OBJECTID:seasonal, geometry)
rm(b4)

saveRDS(b1_filt, paste0(out_fold, "/RISK_32_Polygon_nat_all_BIRD1"))
saveRDS(b2_filt, paste0(out_fold, "/RISK_32_Polygon_nat_all_BIRD2"))
saveRDS(b3_filt, paste0(out_fold, "/RISK_32_Polygon_nat_all_BIRD3"))
saveRDS(b4_filt, paste0(out_fold, "/RISK_32_Polygon_nat_all_BIRD4"))


#### Extract overlapping cells from grid ####

# path for loading polygons
path_poly = "~/predict_vulnerability/Polygons_native_all_bmr/"

# define grids


# world map to filter out oceanic ranges
path_data <- "Z:/THESE/5_Data/Distribution_spatiale/"

land_water1 <- raster(paste0(path_data, "Human_modification/",
                             "gHM_landLakeReservoirOcean300m-0000000000-0000000000.tif"))

land_water2 <- raster(paste0(path_data, "Human_modification/",
                             "gHM_landLakeReservoirOcean300m-0000000000-0000065536.tif"))

land_water3 <- raster(paste0(path_data, "Human_modification/",
                             "gHM_landLakeReservoirOcean300m-0000000000-0000131072.tif"))


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
out_path <- "~/predict_vulnerability/Cells_nat_all_bmr/"

# Computation ----------------

# apply the function to each line of Large spatial polygon df object
# separate per class

#### mammals ####

mam_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_MAM"))
mam_nat_ok <- mam_nat %>%
  dplyr::select(binomial, geometry) %>%
  dplyr::mutate(objectid = 1:nrow(mam_nat))
rm(mam_nat)

mam_IDs <- as.list(mam_nat_ok$objectid)

# for raster at 0.1 degres
df_mam_0.1_mclpply <- pbmclapply(
  mam_IDs,
  function(x){
    poly <- mam_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 10)
saveRDS(df_mam_0.1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_MAM_r01"))
# total SR
mam_df <- bind_rows(df_mam_0.1_mclpply) %>% distinct()
m_all_agg <- mam_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_MAM_r01"))

rm(df_mam_0.1_mclpply, m_all_agg, mam_df, raster0.1, SR_m_as)

# for raster at 1 degres 
df_mam_1_mclpply <- pbmclapply(
  mam_IDs,
  function(x){
    poly <- mam_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 10)
saveRDS(df_mam_1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_MAM_r1"))
# total SR
mam_df <- bind_rows(df_mam_1_mclpply) %>% distinct()
m_all_agg <- mam_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_MAM_r1"))


#### reptiles ####

rept_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_REPT"))
rept_nat_ok <- rept_nat %>%
  dplyr::select(binomial, geometry) %>%
  dplyr::mutate(objectid = 1:nrow(rept_nat))
rm(rept_nat)

rept_IDs <- as.list(rept_nat_ok$objectid)

# for raster at 0.1 degres
df_rept_0.1_mclpply <- pbmclapply(
  rept_IDs,
  function(x){
    poly <- rept_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 10) # Elapsed 59:50
saveRDS(df_rept_0.1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_REPT_r01"))
# total SR
rept_df <- bind_rows(df_rept_0.1_mclpply) %>% distinct()
m_all_agg <- rept_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_REPT_r01"))

rm(df_rept_0.1_mclpply, m_all_agg, rept_df, raster0.1)

# for raster at 1 degres 
df_rept_1_mclpply <- pbmclapply(
  rept_IDs,
  function(x){
    poly <- rept_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 10) # Elapsed 02:59
saveRDS(df_rept_1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_REPT_r1"))
# total SR
rept_df <- bind_rows(df_rept_1_mclpply) %>% distinct()
m_all_agg <- rept_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_REPT_r1"))

#### birds ####

# too heavy for all together
# do it fo each chunk
chunk = 4 # 1, 2, 3 or 4
bird_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD", chunk))

bird_nat_ok <- bird_nat %>%
  dplyr::select(binomial, geometry) %>%
  dplyr::mutate(objectid = 1:nrow(bird_nat))
rm(bird_nat)

bird_IDs <- as.list(bird_nat_ok$objectid)

# for raster at 0.1 degres
df_bird_0.1_mclpply <- pbmclapply(
  bird_IDs,
  function(x){
    poly <- bird_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster0.1)
    return(cells)},
  mc.cores = 5)
saveRDS(df_bird_0.1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_BIRD_r01_chunk", chunk))
# total SR
bird_df <- bind_rows(df_bird_0.1_mclpply) %>% distinct()
m_all_agg <- bird_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_BIRD_r01_chunk", chunk))

rm(df_bird_0.1_mclpply, m_all_agg, bird_df, raster0.1)

# for raster at 1 degres 
df_bird_1_mclpply <- pbmclapply(
  bird_IDs,
  function(x){
    poly <- bird_nat_ok %>%
      dplyr::filter(objectid == x)
    cells <- extract_cells(poly, raster1)
    return(cells)},
  mc.cores = 10)
saveRDS(df_bird_1_mclpply, paste0(out_path, "RISK_32_cells_nat_all_BIRD_r1_chunk", chunk))
# total SR
bird_df <- bind_rows(df_bird_1_mclpply) %>% distinct()
m_all_agg <- bird_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_BIRD_r1_chunk", chunk))


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




SR_m_as<- ggplot(data = m_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_m_as

