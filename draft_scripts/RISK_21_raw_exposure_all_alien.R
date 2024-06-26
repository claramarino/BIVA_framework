# open alien occurrences
# compute first exposure map

rm(list=ls())

library(tidyverse)
library(sp)
library(sf)
library(raster)

# species key
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")


# In a raster, extract alien presence in each cell (for each alien)

#### define extract cells function + raster file ####

# INPUT : a spatial point object and a raster of desired resolution
# OUTPUT : a df containing all the cells with a point, and with nb of pts per cell
extract_cells_pts <- function(pts, rast){
  
  # extract cells with a point inside
  c_pts <- as.data.frame(extract(rast, pts, cellnumbers=TRUE))
  # for each cell, count nb of points inside (density)
  c_pts_dens <- c_pts %>% group_by(cells) %>% summarize(n_pts = n())
  
  # set up species name
  binomial_df <- data.frame(
    new_species = rep(unique(pts$new_species), nrow(c_pts_dens)),
    new_key = rep(unique(pts$new_key), nrow(c_pts_dens)))
  
  # create a df containing:
  df <- bind_cols(
    # xy coordinates from the cells with at least one point inside 
    as.data.frame(xyFromCell(rast, c_pts_dens$cells)),
    # cell name and points density for each cell
    c_pts_dens,
    # name of the species
    binomial_df)
  
  return(df)
}

plot(st_make_grid(what = "centers"), axes = TRUE)
plot(st_make_grid(what = "centers"), axes = TRUE)
plot(st_make_grid(what = "corners"), add = TRUE, col = 'green', pch=3)
sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,0)))))
plot(st_make_grid(sfc, cellsize = .1, square = FALSE))
plot(sfc, add = TRUE)


# Create equal area cell grid


# Using Berhmann's cylindrical equal-area projection
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection




world <- rnaturalearth::countries110
crs(world)
world <- sp::spTransform(world, cea)
# Convert SpatialPolygonsDataframe into sf objects with CEA projection
study_area_sf <- sf::st_transform(sf::st_as_sf(world), crs=cea)

# create grid cells 110 km
grid_110km = st_make_grid(
  x= study_area_sf, cellsize = 110000, what = "polygons", square = F) %>%
  st_sf() 
grid_110km = grid_110km %>%
  mutate(grid_id = 1:nrow(grid_110km)) # add grid ID

# compare projections (make sure both have the same projection)
raster::compareCRS(study_area_sf, grid_110km)

ggplot()+
  geom_sf(data = study_area_sf, fill = "red") +
  geom_sf(data = grid_110km, fill=NA)

# create grid cells 55 km
grid_55km = st_make_grid(
  x= study_area_sf, cellsize = 55000, what = "polygons", square = F) %>%
  st_sf() 
grid_55km = grid_55km %>%
  mutate(grid_id = 1:nrow(grid_55km)) # add grid ID

raster::compareCRS(study_area_sf, grid_110km, grid_55km)

# Transform coordinates of simple feature
points_sf <- sf::st_transform(sf::st_as_sf(points), crs=cea) # certify that data is in the same projection



raster::compareCRS(grid_55km, points_sf)


# count number of points in each grid
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r
grid_sf$n_colli = lengths(st_intersects(honeycomb_grid_sf, test_points))

# remove grid without value of 0 (i.e. no points in side that grid)
honeycomb_count = filter(honeycomb_grid_sf, n_colli > 0)


# initialize raster grid
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


#### Apply to each sp occ file ####

# path for occurrence files
occ_path <- "Output/True_exposure_alien_species/"
occ_files <- list.files(occ_path)

# output folder
out_fold = "Output/Exposure/Exposure_raw/"

# initialise empty list
all_sp_cells_0.1 <- vector(mode = "list", length = length(occ_files))
all_sp_cells_1 <- vector(mode = "list", length = length(occ_files))

for (i in 1:length(occ_files)){
  print(i)
  
  # load alien occurrences
  occ_sp <- readRDS(paste0(occ_path, occ_files[i]))
  
  # extract cells
  df_sp_0.1 <- extract_cells_pts(occ_sp, raster0.1) # resolution = 0.1 deg
  df_sp_1 <- extract_cells_pts(occ_sp, raster1) # resolution = 1 deg
  
  # save cells in list
  all_sp_cells_0.1[[i]] <- df_sp_0.1
  names(all_sp_cells_0.1)[i] <- unique(df_sp_0.1$new_key)
  saveRDS(all_sp_cells_0.1, paste0(out_fold, "RISK_21_raster_cells_308_IAS_r0.1"))
  # save cells in list
  all_sp_cells_1[[i]] <- df_sp_1
  names(all_sp_cells_1)[i] <- unique(df_sp_1$new_key)
  saveRDS(all_sp_cells_1, paste0(out_fold, "RISK_21_raster_cells_308_IAS_r1"))
  
}


#### Explore results ####

# output folder
out_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(out_fold, "RISK_21_raster_cells_308_IAS_r1"))
# combine df for each sp into one df
df_all <- bind_rows(all_sp_cells) 

str(df_all)

# add area correction
area_r1 <- as.data.frame(area(raster1), xy=T) %>% rename(area = layer)
area_r01 <- as.data.frame(area(raster0.1), xy = T) %>% rename(area = layer)

plot(area_r1$area, area_r1$y)

# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())


##### CHECK ALIEN RANGE AND cell number ####

# for mammals and birds, check if consistent with alien range
# from GAVIA and DAMA

# extract area from GAVIA birds
gavia_fold <- "Z:/THESE/5_Data/Distribution_spatiale/Alien_species/GAVIA_rangemaps/"
shp_files <- list.files(gavia_fold)[grepl("*.shp$", list.files(gavia_fold))]

df_area_birds <- data.frame()
for(s in 1:length(shp_files)){
  shp <- st_read(paste0(gavia_fold, shp_files[s]))
  df_area <- data.frame(
    species = shp$Binomial,
    area = st_area(shp)
  )
  df_area_birds <- bind_rows(df_area_birds, df_area)
}
df_area_birds <- df_area_birds %>%
  group_by(species) %>%
  summarise(area = sum(area))

# extract area from DAMA mammals
dama_fold <- "Z:/THESE/5_Data/Distribution_spatiale/Alien_species/DAMA/DataS1/"
shp_files_m <- list.files(dama_fold)[grepl("*.shp$", list.files(dama_fold))]

df_area_mam <- data.frame()
for(s in 1:length(shp_files_m)){
  shp <- st_read(paste0(dama_fold, shp_files_m[s]))
  shp <- st_make_valid(shp)
  df_area <- data.frame(
    species = shp$Binomial,
    area = st_area(shp)
  )
  df_area_mam <- bind_rows(df_area_mam, df_area)
}

df_area_mam <- df_area_mam %>%
  group_by(species) %>%
  summarise(area = sum(area))

# join with tot range
df_area_mb <- bind_rows(
  df_area_birds %>% mutate(class="AVES"), 
  df_area_mam %>% mutate(class = "MAMMALIA")) %>%
  rename(new_species = species)
tot_range_area <- left_join(tot_range, df_area_mb, by = "new_species") %>%
  filter(!is.na(area)) %>%
  mutate(log_cells = log(n_cells_range),
         log_area = log(area/1000000))


library(units)
ggplot(data = tot_range_area, aes(x=log_cells, y=log_area, color = class)) +
  geom_point()+
  geom_smooth(method = lm)

cor.test(tot_range_area$n_cells_range, tot_range_area$area)


##### Compute final raw exposure ####

df_all_range <- left_join(df_all %>% dplyr::select(-new_key), 
                          tot_range,
                           by = "new_species")

# calculate cell related metrics
ias_all_agg <- df_all_range %>%
  group_by(x, y, cells) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_pts), # nb occ per cell (all ias sp)
    nb_occ_med = median(n_pts), # median nb occ per cell
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range)
  ) 

str(ias_all_agg)
summary(ias_all_agg)

# visualize correlations between all metrics
cormat <- round(cor(ias_all_agg),4)
ggcorrplot(cormat)


coast <- rnaturalearth::ne_coastline(returnclass = "sf")

# plot species richness
SR_ias<- ggplot(data = ias_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = coast, alpha = 0.1) +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$new_key))),
       x = "Longitude", y = "Latitude")
SR_ias



