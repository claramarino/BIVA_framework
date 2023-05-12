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
library(pbmcapply)
library(rnaturalearth)

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

#### AMPHIBIANS (for fusias)

# load native ranges from IUCN
amph_all <- st_read(paste0(in_fold, "AMPHIBIANS.shp"))

amph_nat <- amph_all %>% 
  dplyr::filter(origin %in% c(1,2) & presence %in% c(1,2,3,6) & seasonal %in% c(1,2,3)) %>%
  dplyr::select(id_no:seasonal, class:freshwater, geometry) %>%
  dplyr::filter(terrestial == "true" | freshwater == "true") %>%
  dplyr::filter(category %in% c("CR","EN","VU","LC","DD","NT"))


saveRDS(amph_nat, paste0(out_fold, "/RISK_32_Polygon_nat_all_AMPH"))



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

##### Create equal area cell grid for lands ------------

# Using Mollweide equal-area projection as for SDMs
moll <- crs("+proj=moll")
moll

# load worldmap for spatial extent 
worldMap <- ne_countries(scale = "large", type = "countries", returnclass = "sf")
ggplot(worldMap)+geom_sf()
# Convert into sf objects with Mollweide equal-area projection
wp_sf <- sf::st_transform(sf::st_as_sf(worldMap), crs=moll)
ggplot(wp_sf)+geom_sf()
st_bbox(wp_sf)

# create grid cells 110 km
grid_110km = st_make_grid(
  x= wp_sf, cellsize = 110000, what = "polygons", square = F) %>%
  st_sf() 
grid_110km = grid_110km %>%
  mutate(grid_id = 1:nrow(grid_110km)) # add grid ID

# compare projections (make sure both have the same projection)
raster::compareCRS(wp_sf, grid_110km)
# ggplot()+
#   geom_sf(data = wp_sf, fill = "red") +
#   geom_sf(data = grid_110km, fill=NA)

# intersect world and grid for removing non terrestrial cells
inter_110 = st_intersects(wp_sf, grid_110km)
cells110 = unique(c(unlist(inter_110)))
grid_110km_terr <- grid_110km[cells110, ]
# ggplot(grid_110km_terr)+geom_sf()


# create grid cells 55 km
grid_55km = st_make_grid(
  x= wp_sf, cellsize = 55000, what = "polygons", square = F) %>%
  st_sf() 
grid_55km = grid_55km %>%
  mutate(grid_id = 1:nrow(grid_55km)) # add grid ID

raster::compareCRS(grid_110km, grid_55km)
# ggplot()+
#   geom_sf(data = wp_sf, fill = "yellow") +
#   geom_sf(data = grid_55km, fill=NA, alpha=.3)

# intersect world and grid for removing non terrestrial cells
inter_55 = st_intersects(wp_sf, grid_55km)
cells55 = unique(c(unlist(inter_55)))
grid_55km_terr <- grid_55km[cells55, ]
# ggplot(grid_55km_terr)+geom_sf(fill = "red", color = NA)



# save grids for exposure
saveRDS(grid_110km_terr, "Output/RISK_32_grid_110km")
saveRDS(grid_55km_terr, "Output/RISK_32_grid_55km")


##################### Computation of intersection ----------------

# folder for saving outputs
out_path <- "~/predict_vulnerability/Cells_nat_all_bmr/"
# load grids
grid_55km_terr <- readRDS("~/predict_vulnerability/RISK_32_grid_55km")
grid_110km_terr <- readRDS("~/predict_vulnerability/RISK_32_grid_110km")

# apply the function to each line of Large spatial polygon df object
# separate per class

#### mammals ####

mam_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_MAM"))
mam_nat_ok <- st_transform(mam_nat %>%
                             dplyr::select(binomial, geometry) %>%
                             dplyr::mutate(objectid = 1:nrow(mam_nat)),
                           crs = moll) # make sure crs are ok
rm(mam_nat)

mam_IDs <- as.list(mam_nat_ok$objectid)

# test de l'intersection
# test <- st_intersects(mam_nat_ok %>% filter(objectid==mam_IDs[[1]]), 
#                       grid_55km_terr)
# ggplot()+
#   geom_sf(data = mam_nat_ok %>% filter(objectid==mam_IDs[[1]]), fill= "red")+
#   geom_sf(data = grid_55km_terr[unlist(test),], fill = NA)


df_inter_mam_grids <- pbmclapply(
  mam_IDs,
  function(x){
    poly <- mam_nat_ok %>%
      dplyr::filter(objectid == x)
    inter_110 <- st_intersects(poly, grid_110km_terr)
    inter_55 <- st_intersects(poly, grid_55km_terr)
    
    all <- list(cells110 = unique(unlist(inter_110)),
                cells55 = unique(unlist(inter_55)))
    return(all)},
  mc.cores = 10)
# 14 min for mammals

saveRDS(df_inter_mam_grids, paste0(out_path, "RISK_32_cells_nat_all_MAM_110_55"))


rm(mam_nat_ok, df_inter_mam_grids, df_mam_0.1_mclpply, mam_IDs)

#### reptiles ####

rept_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_REPT"))
rept_nat_ok  <- st_transform(rept_nat %>%
                               dplyr::select(binomial, geometry) %>%
                               dplyr::mutate(objectid = 1:nrow(rept_nat)),
                             crs = moll)
rm(rept_nat)

rept_IDs <- as.list(rept_nat_ok$objectid)

# test
# poly <- rept_nat_ok %>%
#   dplyr::filter(objectid == 1)
# inter_110 <- st_intersects(poly, grid_110km_terr)
# cells110 = unique(unlist(inter_110))
# rm(poly, cells110, inter_110)

df_inter_rept_grids <- pbmclapply(
  rept_IDs,
  function(x){
    poly <- rept_nat_ok %>%
      dplyr::filter(objectid == x)
    inter_110 <- st_intersects(poly, grid_110km_terr)
    inter_55 <- st_intersects(poly, grid_55km_terr)
    
    all <- list(cells110 = unique(unlist(inter_110)),
                cells55 = unique(unlist(inter_55)))
    return(all)},
  mc.cores = 10)

saveRDS(df_inter_rept_grids, paste0(out_path, "RISK_32_cells_nat_all_REPT_110_55"))

rm(rept_nat_ok, df_inter_rept_grids, df_rept_0.1_mclpply, rept_IDs)



#### amphibians ####

amph_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_AMPH"))
amph_nat_ok  <- st_transform(amph_nat %>%
                               dplyr::select(binomial, geometry) %>%
                               dplyr::mutate(objectid = 1:nrow(amph_nat)),
                             crs = moll)
rm(amph_nat)

amph_IDs <- as.list(amph_nat_ok$objectid)

# test
# poly <- amph_nat_ok %>%
#   dplyr::filter(objectid == 1)
# inter_110 <- st_intersects(poly, grid_110km_terr)
# cells110 = unique(unlist(inter_110))
# rm(poly, cells110, inter_110)

df_inter_amph_grids <- pbmclapply(
  amph_IDs,
  function(x){
    poly <- amph_nat_ok %>%
      dplyr::filter(objectid == x)
    inter_110 <- st_intersects(poly, grid_110km_terr)
    inter_55 <- st_intersects(poly, grid_55km_terr)
    
    all <- list(cells110 = unique(unlist(inter_110)),
                cells55 = unique(unlist(inter_55)))
    return(all)},
  mc.cores = 10)

saveRDS(df_inter_amph_grids, paste0(out_path, "RISK_32_cells_nat_all_AMPH_110_55"))


#### birds ####

# too heavy for all together
# do it fo each chunk
chunk = 1 # 1, 2, 3 or 4
bird_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD", chunk))

bird_nat_ok  <- st_transform(bird_nat %>%
                               dplyr::select(binomial, geometry) %>%
                               dplyr::mutate(objectid = 1:nrow(bird_nat)),
                             crs = moll)
rm(bird_nat)

bird_IDs <- as.list(bird_nat_ok$objectid)


# test
poly <- bird_nat_ok %>%
  dplyr::filter(objectid == 1)
inter_110 <- st_intersects(poly, grid_110km_terr)
cells110 = unique(unlist(inter_110))


df_inter_bird_grids <- pbmclapply(
  bird_IDs,
  function(x){
    poly <- bird_nat_ok %>%
      dplyr::filter(objectid == x)
    inter_110 <- st_intersects(poly, grid_110km_terr)
    inter_55 <- st_intersects(poly, grid_55km_terr)
    
    all <- list(cells110 = unique(unlist(inter_110)),
                cells55 = unique(unlist(inter_55)))
    return(all)},
  mc.cores = 10)

saveRDS(df_inter_bird_grids,
        paste0(out_path, "RISK_32_cells_nat_all_BIRD_110_55_chunk", chunk))




#############################


world <- rnaturalearth::countries110

study_area_sf <- sf::st_transform(sf::st_as_sf(world), crs=moll)
ggplot(study_area_sf)+geom_sf()


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



# total SR
mam_df <- bind_rows(df_mam_0.1_mclpply) %>% distinct()
m_all_agg <- mam_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_MAM_r01"))


# total SR
rept_df <- bind_rows(df_rept_0.1_mclpply) %>% distinct()
m_all_agg <- rept_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_REPT_r01"))

# total SR
bird_df <- bind_rows(df_bird_0.1_mclpply) %>% distinct()
m_all_agg <- bird_df %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())
saveRDS(m_all_agg, paste0(out_path, "RISK_32_SR_tot_BIRD_r01_chunk", chunk))