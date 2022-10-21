# Explore GBIF coverage data
# from meyer et al 2015


rm(list=ls())

library(tidyverse)
library(sp)
library(sf)
library(ggthemes)

path <- "Z:/THESE/5_Data/Distribution_spatiale/GBIF_data_bias_Meyer/"
list.files(path)
amph <- read.csv2(paste0(path, "Meyer_etal_data_amph_110km.csv"))
mam <- read.csv2(paste0(path, "Meyer_etal_data_mamm_110km.csv"))
bird <- read.csv2(paste0(path, "Meyer_etal_data_aves_110km.csv"))

str(amph)

ggplot(amph) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = InventoryCompleteness)) +
  scale_fill_gradient(
    low = "red", 
    high = "green")+
  coord_map("moll")+
  theme_map()


ggplot(mam) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = InventoryCompleteness)) +
  scale_fill_gradient(
    low = "red", 
    high = "green")+
  theme_classic()

ggplot(bird) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = InventoryCompleteness)) +
  scale_fill_gradient(
    low = "red", 
    high = "green")+
  coord_map("moll")+
  theme_map()



all <- full_join(
  full_join(bird %>% dplyr::select(Longitude, Latitude, InventoryCompleteness) %>%
              rename(x=Longitude, y = Latitude, Compl_bird = InventoryCompleteness),
            mam %>% dplyr::select(Longitude, Latitude, InventoryCompleteness) %>%
              rename(x=Longitude, y = Latitude, Compl_mam = InventoryCompleteness)),
  amph %>% dplyr::select(Longitude, Latitude, InventoryCompleteness) %>%
    rename(x=Longitude, y = Latitude, Compl_amph = InventoryCompleteness)) %>%
  # calculate mean completeness across the 3 taxa 
  mutate(Compl_mean = rowMeans(bind_cols(Compl_bird, Compl_mam, Compl_amph), na.rm = T))


plot(all$Compl_mean, all$Compl_mam)
hist(all$Compl_mean)

hist(all$Compl_mam)
hist(all$Compl_amph)
hist(all$Compl_bird)


# attribute a completeness value to each grid of our rasters
library(raster)

all_coord <- all %>%
  st_as_sf(
    coords = c("x", "y"),
    agr = "constant",
    crs = "+proj=longlat +datum=WGS84",
    stringsAsFactors = FALSE,
    remove = TRUE)

# Create raster
raster1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                  resolution = 1, vals=NULL)

# rasterize points
compl_rast <- raster::rasterize(all_coord, raster1, "Compl_mean", fun = mean)

plot(compl_rast)

# pbm of latitudes without info (due to equal area proj)
# correct cells with vertical neighbors
w <- matrix(c(0,0,0,1,0,1,0,0,0), nrow =3)
w
compl_rast_foc <- focal(compl_rast, w, max, na.rm=T, NAonly=TRUE, pad=F)
plot(compl_rast_foc)

