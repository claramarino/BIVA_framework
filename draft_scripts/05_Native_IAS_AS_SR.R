rm(list=ls())

library(tidyverse)
library(rredlist)
library(raster)
library(rnaturalearth)
library(sp)
library(rgeos)
library(maptools)
library(sf)
library(rgdal)
library(RColorBrewer)

# Load species threats
setwd("Z:/THESE/6_Projects/predict_vulnerability")
df_threat <- readRDS("Output/Data_clean/03_df_all_threat_assoc_bmr")

# list_df_threat_ft <- readRDS("Output/Data_clean/03_df_all_threat_assoc_bmr_with_traits")
# 
# # extract species list for each class
# list_sp_as_2more <- lapply(list_df_threat_ft, function(x){
#   x %>% filter(Assoc_ias_2more == 1) %>% pull(binomial)})
# lapply(list_df_threat_ft, nrow)
# 10310 + 5502 + 4901
# # when taking only species with traits, 3200 species are missing from df
# # I prefer taking all the species 

sp_list <- df_threat %>% 
  filter(Assoc_ias_2more==1)%>%
  pull(binomial)

#--------------------------- Range size ----------------------------

#------------------- Birds

# my_path <- "Z:/THESE/5_Data/Distribution_spatiale"
# setwd(paste(my_path,"/1_MAJ_BIRDS_Birdlife", sep=""))
# # Load the 4 shapefiles
# birds_all_1 <- readOGR("Birds_A_2600_repare.shp")
# birds_all_2 <- readOGR("Birds_A_3400_repare.shp")
# birds_all_3 <- readOGR("Birds_B_5800_repare.shp")
# birds_all_4 <- readOGR("Birds_C_5766_repare.shp")
# 
# birds_in_list1 <- subset(birds_all_1, binomial %in% sp_list)
# birds_in_list2 <- subset(birds_all_2, binomial %in% sp_list)
# birds_in_list3 <- subset(birds_all_3, binomial %in% sp_list)
# birds_in_list4 <- subset(birds_all_4, binomial %in% sp_list)
# 
# birds_in_list_all <- rbind(rbind(birds_in_list1, birds_in_list2),
#                    rbind(birds_in_list3, birds_in_list4))
# # filter for native range only
# birds_in_list_all_native <- subset(birds_in_list_all, 
#                                    presence==1 & (origin==1 | origin==2))
# 
# setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Polygons_IAS_AS")
# saveRDS(birds_in_list_all_native, "05_polygons_birds_IAS_as_2more")

# # do all sp have a polygon?
# sp_poly <- unique(birds_ias_as$binomial)
# table(df_threat %>% filter(Assoc_ias_2more == 1) %>% pull(Class))
# # 1 sp is missing
# setdiff(df_threat %>% filter(Assoc_ias_2more == 1 & Class=="Aves") %>% pull(binomial),
#         sp_poly)
# Grus Antigone... => à refaire tourner en changeant son nom au début

# # Load oceans from natural earth
# oceans <- readOGR("Z:/THESE/5_Data/Distribution_spatiale/Shp_oceans_natural_earth_10m/ne_10m_ocean.shp")
# plot(oceans, col="blue")
# 
# # compute intersection between birds and oceans (takes ~15h)
# b_oc <- raster::intersect(birds_ias_as, oceans)
# saveRDS(b_oc, "05_intersect_ocean_birds")

#------------------- Mammals

# my_path <- "Z:/THESE/5_Data/Distribution_spatiale/"
# mam_all <- readOGR(paste0(my_path,"0_MAJ_IUCN_2020/Mam_repare.shp"))
# 
# # select all mammals, that are native & extant
# mam_in_list_nat <- subset(mam_all, binomial %in% sp_list &
#                         presence==1 & (origin==1 | origin==2))
# 
# setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Polygons_IAS_AS")
# saveRDS(mam_in_list_nat, "05_polygons_mam_IAS_as_2more")

# # do all sp have a polygon?
# sp_poly <- unique(mam_in_list_nat$binomial)
# table(df_threat %>% filter(Assoc_ias_2more == 1) %>% pull(Class))
# # yes

#------------------- Reptiles

# rept_all <- readOGR(
#   paste0(my_path,"REPTILE_Distrib_GARD.1_dissolved_ranges/Reptiles_GARD_repare.shp"))
# rept_in_list <- subset(rept_all, Binomial %in% sp_list)
# setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Polygons_IAS_AS")
# saveRDS(rept_in_list, "05_polygons_rept_IAS_as_2more")

# # do all sp have a polygon?
# rept_poly <- unique(rept_in_list$Binomial)
# # no 20 species have no polygons

# ----------------------- Load final polygons all classes

# setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Polygons_IAS_AS")
# 
# birds_ias_as <- readRDS("05_polygons_birds_IAS_as_2more")
# mam_ias_as <- readRDS("05_polygons_mam_IAS_as_2more")
# rept_ias_as <- readRDS("05_polygons_rept_IAS_as_2more")
# birds_oceans <- readRDS("05_intersect_ocean_birds")


#---------------------------- Rasterize richness ----------------

# # define raster
# raster1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, 
#                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
#                   resolution = 1, vals=NULL)
# raster0.5 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, 
#                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", 
#                   resolution = 0.5, vals=NULL)
# 
# # folder for saving outputs
# setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Raster_IAS_AS")


# ------------ for birds

# # takes 1 h 20 min 
# Sys.time()
# r.sp.b.0.5 <- rasterize(x = birds_ias_as, 
#                   y = raster0.5, 
#                   field = "binomial",
#                   fun = function (x, ...) length(unique(na.omit(x)))) 
# # pbm : on a tous les oiseaux marins, donc leur aire de distrib est dans les océans...
# saveRDS(r.sp.b.0.5, "05_raster_richness_birds_ias_as_0_5")
# 
# Sys.time()
# # rasterize birds from ocean
# r.sp.boc.0.5 <- rasterize(x = birds_oceans, 
#                         y = raster0.5, 
#                         field = "binomial",
#                         fun = function (x, ...) length(unique(na.omit(x)))) 
# saveRDS(r.sp.boc.0.5, "05_raster_richness_birds_ias_as_ocean_0_5")
# Sys.time()

# # compute species richness difference between terrestrial and oceans
# # change NA to zéro in ocean raster
# r.sp.boc.0.5[is.na(r.sp.boc.0.5[])] <- 0 
# # for large rasters, don't perform raster math but use overlay function
# # here okay with raster math
# r.sp.bland.0.5 <- r.sp.b.0.5 - r.sp.boc.0.5
# plot(r.sp.bland.0.5)
# saveRDS(r.sp.bland.0.5, "05_raster_richness_birds_ias_as_land_0_5")

# ------------ for mammals

# takes ~1 min for resolution =1
# takes ~2 min for resolution =0.5

# r.sp.m.0.5 <- rasterize(x = mam_ias_as, 
#                   y = raster0.5, 
#                   field = "binomial",
#                   fun = function (x, ...) length(unique(na.omit(x)))) 
# saveRDS(r.sp.m.0.5, "05_raster_richness_mam_ias_as_0_5")
# 
# #plot(r.sp.m)
# Sys.time()

# # ------------ for reptiles

# # takes ~30 sec 
# r.sp.r.0.5 <- rasterize(x = rept_ias_as, 
#                   y = raster0.5, 
#                   field = "Binomial",
#                   fun = function (x, ...) length(unique(na.omit(x)))) 
# #plot(r.sp.r.0.5)
# 
# saveRDS(r.sp.r.0.5, "05_raster_richness_rept_ias_as_0_5")
# 
# Sys.time()


#####################################################
# Plot species richness
#####################################################

# load rasters
setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Raster_IAS_AS")
birds05 <- readRDS("05_raster_richness_birds_ias_as_land_0_5")
mam05 <- readRDS("05_raster_richness_mam_ias_as_0_5")
rept05 <- readRDS("05_raster_richness_rept_ias_as_0_5")

# world map
worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")


# str(birds05)
# str(mam05)
# str(rept05)


# convert raster to df
mam05_df <-
  as.data.frame(mam05, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() %>%
  #--- change the variable names ---#
  rename(SR = layer) %>%
  mutate(Class = "Mammalia")

rept05_df <-
  as.data.frame(rept05, xy = TRUE) %>%
  na.omit() %>%
  rename(SR = layer) %>%
  mutate(Class = "Reptilia")

birds05_df <-
  as.data.frame(birds05, xy = TRUE) %>%
  na.omit() %>%
  rename(SR = layer) %>%
  mutate(Class = "Aves") %>%
  filter(SR>0) # attention, voir d'où vient le pbm

head(mam05_df)
head(rept05_df)
head(birds05_df)

hist(mam05_df$SR)
ggplot(mam05_df) +
  geom_point(aes(x =SR, y=y), position = "jitter")

hist(birds05_df$SR)
ggplot(birds05_df) +
  geom_point(aes(x =SR, y=y), position = "jitter")

hist(rept05_df$SR)
ggplot(rept05_df) +
  geom_point(aes(x =SR, y=y), position = "jitter")




SR_mam <- ggplot(data = mam05_df) +
  geom_raster(aes(x = x, y = y, fill = SR)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Taxonomic richness of mammals associated to one of the target IAS",
       subtitle = paste0("Total number of species: ",
                         nrow(df_threat %>% filter(Class == "Mammalia" & Assoc_ias_2more == 1))),
       x = "Longitude", y = "Latitude") 
SR_mam


SR_rept <- ggplot(data = rept05_df) +
  geom_raster(aes(x = x, y = y, fill = SR)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Taxonomic richness of reptiles associated to one of the target IAS",
       subtitle = paste0("Total number of species: ",
                         nrow(df_threat %>% filter(Class == "Reptilia" & Assoc_ias_2more == 1))),
                         x = "Longitude", y = "Latitude")

SR_rept

SR_birds <- ggplot(data = birds05_df) +
  geom_raster(aes(x = x, y = y, fill = SR)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Taxonomic richness of birds associated to one of the target IAS",
       subtitle = paste0("Total number of species: ",
                         nrow(df_threat %>% filter(Class == "Aves" & Assoc_ias_2more == 1))),
       x = "Longitude", y = "Latitude")

SR_birds



bmr_as <- bind_rows(mam05_df, rept05_df, birds05_df)

str(bmr_as)

bmr_as_agg <- bmr_as %>%
  group_by(x, y) %>%
  summarise(SR_tot = sum(SR))


SR_bmr_as<- ggplot(data = bmr_as_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Taxonomic richness of tetrapods associated to one of the target IAS",
       subtitle = paste0("Total number of species (birds, mammals, reptiles): ",
                         nrow(df_threat %>% filter(Assoc_ias_2more == 1))),
       x = "Longitude", y = "Latitude")

SR_bmr_as





### brouillon

ggplot(data = worldMap) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  ggtitle("World map", 
          subtitle = paste0("(", length(unique(worldMap$admin)), " countries)"))



# grid
wgrid <- st_make_grid(worldMap, cellsize = 0.5)
class(wgrid)

wgrid_large <- st_make_grid(worldMap, cellsize = 2)
ggplot(data = wgrid_large) +
  geom_sf() 

st_crs(wgrid)





plot(r.sp)
