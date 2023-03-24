# open alien occurrences
# compute first exposure map

rm(list=ls())

library(tidyverse)
library(sp)
library(sf)
library(raster)

# species key
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")


# In an equal-area grid, extract alien presence in each cell (for each alien)

#### define grids + function for extracting cell points ####

##### Create equal area cell grid

# Using Berhmann's cylindrical equal-area projection (CEA)
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load worldmap for spatial extent 
world <- rnaturalearth::countries110
# Convert into sf objects with CEA projection
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

range(grid_55km$grid_id)

raster::compareCRS(grid_110km, grid_55km)

##### Define function for extraction

# INPUT : a spatial point object and a grid of desired resolution
# OUTPUT : a df containing all the cells with a point, and with nb of pts per cell
extract_cells_pts <- function(pts, grid){

  # attribute the good crs 
  pts_cea<- sf::st_transform(sf::st_as_sf(pts), crs=cea)
  
  if (raster::compareCRS(pts_cea, grid)){
    # count nb of points per grid cell
    grid$n_occ <- lengths(st_intersects(grid, pts_cea))
    # keep only cells with occurrences
    grid_fin <- grid %>% filter(n_occ>0) %>%
      mutate(new_species = unique(pts_cea$new_species),
             new_key = unique(pts_cea$new_key)) %>%
      st_drop_geometry()    
  } else {print("Pbm crs")} 

}



#### Apply to each sp occ file ####

# path for occurrence files
occ_path <- "Output/True_exposure_alien_species/"
occ_files <- list.files(occ_path)

# output folder
out_fold = "Output/Exposure/Exposure_raw/"

# initialise empty df
all_sp_cells_55 <- all_sp_cells_110 <- all_sp_cells_110_nofilter <- data.frame()


for (i in 1:length(occ_files)){
  print(i)
  
  # load alien occurrences
  occ_sp <- readRDS(paste0(occ_path, occ_files[i]))
  
  occ_sp_55 <- occ_sp %>% dplyr::filter(coordinateUncertaintyInMeters < 22500)
  # occ_sp_110 <- occ_sp %>% dplyr::filter(coordinateUncertaintyInMeters < 55000)
  
  # extract cells
  df_sp_55 <- extract_cells_pts(occ_sp_55, grid_55km) # resolution = 55km x 55km (0.5 deg)
  # df_sp_110 <- extract_cells_pts(occ_sp_110, grid_110km) # resolution = 110km x 110km (1 deg)
  # df_sp_110_nofilter <- extract_cells_pts(occ_sp, grid_110km) # resolution = 110km x 110km (1 deg)
  # 
  # save cells in dataframes
  all_sp_cells_55 <- bind_rows(all_sp_cells_55, df_sp_55)
  saveRDS(all_sp_cells_55, paste0(out_fold, "RISK_21_grid_cells_310_IAS_55km"))
  # # save cells in list
  # all_sp_cells_110 <- bind_rows(all_sp_cells_110, df_sp_110)
  # saveRDS(all_sp_cells_110, paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km"))
  # # save cells in list
  # all_sp_cells_110_nofilter <- bind_rows(all_sp_cells_110_nofilter, df_sp_110_nofilter)
  # saveRDS(all_sp_cells_110_nofilter, paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km_nofilter"))
  
}
st_bbox(grid_110km)
st_bbox(grid_55km)

#### Explore results ####

# output folder
out_fold = "Output/Exposure/Exposure_raw/"

# df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km_nofilter"))
# df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_110km"))
df_all <- readRDS(paste0(out_fold, "RISK_21_grid_cells_310_IAS_55km"))

range(grid_55km$grid_id)
range(df_all$grid_id)


# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())

hist(tot_range$n_cells_range)

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
mod <- lm(n_cells_range ~ area + class, data = tot_range_area)
summary(mod)




##### Compute final raw exposure ####

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

str(ias_all_agg)
summary(ias_all_agg)

# visualize correlations between all metrics
cormat <- round(cor(ias_all_agg),4)
library(ggcorrplot)
ggcorrplot(cormat)

# plot alien species richness 

# Using Berhmann's cylindrical equal-area projection (CEA)
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load worldmap for spatial extent 
world <- rnaturalearth::countries110
# Convert into sf objects with CEA projection
study_area_sf <- sf::st_transform(sf::st_as_sf(world), crs=cea)
# create grid cells 110 km
grid_110km = st_make_grid(
  x= study_area_sf, cellsize = 110000, what = "polygons", square = F) %>%
  st_sf() 
grid_110km = grid_110km %>%
  mutate(grid_id = 1:nrow(grid_110km))

grid_55km = st_make_grid(
  x= study_area_sf, cellsize = 55000, what = "polygons", square = F) %>%
  st_sf() 
grid_55km = grid_55km %>%
  mutate(grid_id = 1:nrow(grid_55km)) # add grid ID

range(grid_55km$grid_id)

ias_agg_grid110 <- st_as_sf(left_join(ias_all_agg, grid_110km))
ias_agg_grid55 <- st_as_sf(left_join(ias_all_agg, grid_55km))


# plot species richness
SR_ias<- ggplot(data = ias_agg_grid55) +
  geom_sf(aes(fill = SR_tot), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1) +
  geom_sf(data = study_area_sf, fill = NA, color = "grey70") +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$new_key))),
       x = "Longitude", y = "Latitude")
SR_ias


ggplot(grid %>% filter(n_occ>0), aes(fill=n_occ))+
  geom_sf(color = NA)+
  scale_fill_gradient(
    low = "white", 
    high = "red")



########################
# compare with previous method (no equal carea + no filter on coord uncertainty)


# output folder
out_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(out_fold, "RISK_21_raster_cells_308_IAS_r1"))
# combine df for each sp into one df
df_all2 <- bind_rows(all_sp_cells) 

# for each ias, calculate max range
tot_range2<- df_all2 %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())

range_all <- left_join(tot_range, tot_range2, by ="new_species")
ggplot(range_all, aes(x= n_cells_range.x, y = n_cells_range.y))+
  geom_point() + geom_abline(slope = 1, intercept = 0)
