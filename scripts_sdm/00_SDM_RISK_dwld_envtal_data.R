# environmental data for SDMs projected exposure
rm(list = ls())
library(readxl)
library(downloader)
library(terra)
library(raster)

# file path
sdm_path <- "Z:/R/sdmrisk/"

##### CHELSA download #####

for(bioclim in sprintf("%02d", 1:19)){
  # download baseline climatic data
  download(
    paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_", bioclim, ".tif"), 
    paste0(sdm_path, "data/raw/bioclim/CHELSA_baseline_bio", bioclim, ".tif"),
    mode = "wb")
}

##### Aggregate LU to decrease the resolution #####

# open raw data from copernicus (100m x 100m)
lu_cate <- terra::rast(paste0(
  sdm_path, 
  "data/raw/lu/",
  "PROBAV_LC100_global_v3.0.1_2015-base_Discrete-Classification-map_EPSG-4326.tif"))

object.size(lu_cate)
lu_cate
ext(lu_cate)

# define lu_cate code following Copernicus metadata
lu_code <- data.frame(
  orig_code = c(0,111,113,112,114,115,116,121,123,122,124,125,126,
                20,30,90,100,60,40,50,70,80,200),
  new_code = c(0,110,110,110,110,110,110,110,110,110,110,110,110,
               20,30,80,100,60,40,50,70,80,200),
  landcover = c("NA","forest","forest","forest","forest","forest","forest",
                "forest","forest","forest","forest","forest","forest",
                "shrubs","grass", "water","moss_lichen","bare","crop",
                "built","snow_ice","water", "sea")
)

# get unique codes with new categories
codes <- as.character(unique(lu_code$new_code))
aggreg = 27 # define the aggregation factor (multiple of 9, 7, 8)

# x in 1:20, y in 1:16
k=0
for (x in 1:20){
  print(Sys.time())
  for(y in 1:16){
    k=k+1
    # crop extent for 320 squares on the world
    crp <- terra::crop(lu_cate, 
                       ext(-180 + 18*(x-1),
                           -180 + 18*x, 
                           80 - 9*y, 
                           80 - 9*(y-1)))
    # change values with adequate categories
    crp_sub <- subst(crp, lu_code$orig_code, lu_code$new_code)
    # segregate by class => each land cover is a dummy layer with 0/1 value
    slu_crp <- segregate(crp_sub)
    # take the mean across 27*27 cells for having the proportion of each land use
    prop <- aggregate(slu_crp, aggreg, mean)
    
    # create an empty raster
    # for stacking all layers of land cover in the same order for each chunck
    # necessary step for the final merging
    all_lay <- prop[[1]]
    names(all_lay) <- "base"
    all_lay[] <- NA
    all_lay
    
    for (i in codes){
      if (i %in% names(prop)){
        to_bind <- prop[[i]]
      } else {
        to_bind <- prop[[1]]
        names(to_bind) <- i
        to_bind[] <- 0
      }
      all_lay <- c(all_lay, to_bind)
    }
    
    # save aggregated raster 
    saveRDS(all_lay, paste0(sdm_path, "data/output/lu/chunks_agg/", 
                            "lu_prop_agg_", aggreg, "_chunck_", x, "_", y))
    
    print(paste0(k, "/320 done"))
  }
  print(Sys.time())
}


# open all chunks in a list
chunks <- list.files(paste0(sdm_path, "data/output/lu/chunks_agg/"), full.names = T)
rlist <- list()
for (i in 1:length(chunks)){
  rlist <- c(rlist, list(readRDS(chunks[i]) ))
}

# create a spatial collection
rsrc <- sprc(rlist)
# merge them all into one single raster 
m <- mosaic(rsrc)
names(m) <- c("base", unique(lu_code$landcover))
#plot(m)

m_ok <- m[[c("NA","forest","shrubs","grass","water","moss_lichen","bare","crop","built","snow_ice")]]
# replace cells that have NA with true NAs
m_ok[m_ok[["NA"]]==1] <- NA
#plot(m_ok[["built"]])

m_ok[["NA"]] <- NULL
m_ok

file_names <- paste0(sdm_path, "data/output/lu/lu_prop_agg_", 
                     aggreg, "_mosaic_", names(m_ok), ".tif")
writeRaster(m_ok, file_names, overwrite=TRUE)

# save final raster
saveRDS(m_ok, paste0(sdm_path, "data/output/lu/", 
                        "lu_prop_agg_", aggreg, "_mosaic.rds"))


##### Target group for correcting sampling bias #####


# folder
fold = "Output/Exposure/Exposure_raw/"
df_all <- readRDS(paste0(fold, "RISK_21_grid_cells_310_IAS_55km"))
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")
gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")

library(tidyverse)

info <- left_join(sp_info, gbif_taxo_2) %>% 
  dplyr::select(new_key:freshwater, kingdom:order,class) %>%
  distinct()

table(info$class)
table(info$kingdom)

table(info$phylum)

# 23 arthropodes + 3 molluscs + 1 ctenophore
# 191 vertebrates dont 17 poissons
# 110 plants

# calculate richness for each group
df_all_info <- left_join(
  df_all %>% dplyr::select(-new_key),
  info %>% mutate(group = case_when(
    kingdom=="Plantae" ~ "plant",
    phylum %in% c("Arthropoda", "Mollusca","Ctenophora") ~ "invert",
    phylum == "Chordata" & class == "Actinopterygii" ~ "vert_fish",
    phylum == "Chordata" & class != "Actinopterygii" ~ "vert_terr"
    )),
  by = "new_species")

# calculate cell related metrics
SR_all_agg <- df_all_info %>%
  group_by(grid_id) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_occ), # nb occ per cell (all ias sp)
  ) 

SR_group_agg <- df_all_info %>%
  group_by(grid_id, group) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_occ), # nb occ per cell (all ias sp)
  )

SR_all_lay <- left_join(SR_all_agg, SR_group_agg %>%
                          pivot_wider(names_from = group, 
                                      values_from = SR_tot:nb_occ_tot))

# get coordinates of grid_id 

# Using Berhmann's cylindrical equal-area projection (CEA)
cea<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load worldmap for spatial extent 
library(sf)
world <- rnaturalearth::countries110
# Convert into sf objects with CEA projection
study_area_sf <- sf::st_transform(sf::st_as_sf(world), crs=cea)
st_bbox(study_area_sf)

# create grid cells 55 km
grid_55km = st_make_grid(
  x= study_area_sf, cellsize = 55000, what = "polygons", square = F) %>%
  st_sf() 
grid_55km = grid_55km %>%
  mutate(grid_id = 1:nrow(grid_55km))
range(grid_55km$grid_id)
range(SR_all_lay$grid_id)



# check the differences in SR between groups for contribution to total SR
ggplot(left_join(SR_group_agg %>% rename(SR = SR_tot), 
                 SR_all_agg %>% select(grid_id:SR_tot))) + 
  geom_point(aes(x=SR_tot, y=SR, color = group))+
  geom_smooth(aes(x=SR_tot, y=SR, color = group), method = "lm")+
  ylab("SR by taxonomic group")
  

# sounds like we can make 3 groups:
# vertebrates_terr, plants, and [invert-vert_fish]

grid_sr <- left_join(grid_55km, SR_all_lay) %>%
  # carreful to sum even if there is NA in one column 
  mutate(SR_fish_invert = replace_na(SR_tot_vert_fish, 0) + replace_na(SR_tot_invert, 0))%>%
  filter(!is.na(SR_tot))
# replace 0 values with NAs in the column SR_fish_invert
grid_sr$SR_fish_invert[grid_sr$SR_fish_invert==0] <- NA


colnames(grid_sr)

class(grid_sr)

# richness all exotic species 
ggplot(grid_sr)+
  geom_sf(aes(fill = SR_tot), color=NA)+
  geom_sf(data = study_area_sf, fill = NA, color = "grey70") +
  scale_fill_viridis_c(option="magma", direction = -1) +
  theme_classic()

# richness exotic terrestrial vertebrates 
ggplot(grid_sr %>% filter(!is.na(SR_tot_vert_terr)))+
  geom_sf(data = study_area_sf, fill = NA, color = "grey70") +
  geom_sf(aes(fill = SR_tot_vert_terr), color=NA)+
  scale_fill_viridis_c(option="magma", direction = -1) +
  theme_classic()
# richness exotic plants
ggplot(grid_sr %>% filter(!is.na(SR_tot_plant)))+
  geom_sf(aes(fill = SR_tot_plant), color=NA)+
  scale_fill_viridis_c(option="magma", direction = -1) +
  geom_sf(data = study_area_sf, fill = NA, color = "grey70") +
  theme_classic()

# richness exotic invertebrates and fish
ggplot(grid_sr %>% filter(!is.na(SR_fish_invert)))+
  geom_sf(aes(fill = SR_fish_invert), color=NA)+
  scale_fill_viridis_c(option="magma", direction = -1) +
  geom_sf(data = study_area_sf, fill = NA, color = "grey70") +
  theme_classic()


# plants vs vert
ggplot(grid_sr)+
  geom_point(aes(x=SR_tot_vert_terr, y=SR_tot_plant))

# plants vs invert
ggplot(grid_sr)+
  geom_point(aes(x=SR_fish_invert, y=SR_tot_plant))

# vert vs invert
ggplot(grid_sr)+
  geom_point(aes(x=SR_tot_vert_terr, y=SR_fish_invert))

range(na.omit(grid_sr$SR_tot_vert_terr))
range(na.omit(grid_sr$SR_tot_plant))
range(na.omit(grid_sr$SR_fish_invert))

# create a raster
# containing sp richness of all exotic species in the db
# and one layer per big class
# plants, terr vert, freshwater vert, invert, plants

# same resolution as climate
chelsa <- rast("Z:/R/sdmrisk/data/raw/bioclim/CHELSA_baseline_bio01.tif")
chelsa_10km <- terra::aggregate(chelsa, fact = 10)
chelsa_10km

# transform polygons in the adequate crs
grid_sr_84 <- sf::st_transform(grid_sr, crs=crs(chelsa_10km)) %>%
  mutate(area = st_area(geometry))
st_bbox(grid_sr)
st_bbox(grid_sr_84)
st_geometry(grid_sr_84)

ggplot(grid_sr_84)+
  geom_sf(aes(fill = as.numeric(area)), color=NA)+
  ggtitle("le pont du monde")

sr_v <- normalize.longitude(vect(grid_sr_84))
crs(sr_v)
crs(chelsa_10km)

r_c <- c(
  terra::rasterize(sr_v, chelsa_10km, field = "SR_tot", fun="sum"),
  terra::rasterize(sr_v, chelsa_10km, field = "SR_tot_vert_terr", fun="sum"),
  terra::rasterize(sr_v, chelsa_10km, field = "SR_tot_plant", fun="sum"),
  terra::rasterize(sr_v, chelsa_10km, field = "SR_fish_invert", fun="sum"))

r_c
plot(r_c)









############################## BROUILLON ##############################"

lu <- data.frame(
  ID = c(0,111,113,112,114,115,116,121,123,122,124,125,126,20,30,90,100,60,40,50,70,80,200),
  landcover = c("NA", "forest_closed","forest_closed","forest_closed","forest_closed",
                "forest_closed", "forest_closed", "forest_open", "forest_open", 
                "forest_open", "forest_open", "forest_open", "forest_open",
                "shrubs","grass", "wetland","moss_lichen","bare","crop",
                "built","snow_ice","water_perm", "sea"), 
  code = c("NA","forest","forest","forest","forest","forest","forest",
           "forest","forest","forest","forest","forest","forest",
           "shrubs","grass", "water","moss_lichen","bare","crop",
           "built","snow_ice","water", "sea"))


levels(lu_cate1) <- lu
lu_cate1
# save forest 2.5
plot(forest)

forest

sum(is.na(forest))

cellStats(forest)


#lu_catef1 <- ratify(lu_cate1)

levels(lu_cate1) <- lu
is.factor(lu_cate1)
values(lu_cate1)[1:100]

levels(lu_cate1)
attributes(lu_cate1)



r <- round(r)
f <- as.factor(r)
is.factor(f)
summary(f)

names(lu_cate1) 
# aggregate 

forest_1 <- lu_cate1
forest_1[forest_1[] == "" ] = NA


r <- raster(nrow=10, ncol=10)
values(r) <- runif(ncell(r)) * 10
is.factor(r)

r <- round(r)
f <- as.factor(r)
is.factor(f)

x <- levels(f)[[1]]
x
x$code <- letters[10:20]
levels(f) <- x
levels(f)
f

r <- raster(nrow=10, ncol=10)
values(r) = 1
r[51:100] = 2
r[3:6, 1:5] = 3
r <- ratify(r)

rat <- levels(r)[[1]]
rat$landcover <- c("Pine", "Oak", "Meadow")
rat$code <- c(12,25,30)
levels(r) <- rat
r

# extract values for some cells
i <- extract(r, c(1,2, 25,100))
i
# get the attribute values for these cells
factorValues(r, i)

# write to file:
# rr <- writeRaster(r, rasterTmpFile(), overwrite=TRUE)
# rr

# create a single-layer factor 
x <- deratify(r, "landcover")
x
is.factor(x)
levels(x)






##### Other datasets #####

# Aggregate chelsa to 10km to rasterize ports & airports
chelsa <- raster("Z:/R/sdmrisk/data/raw/bioclim/CHELSA_baseline_bio01.tif")
chelsa_10km <- raster::aggregate(chelsa, 
                                 fact = 10)

# Roads (motorways)
library(sf)
roads <- st_read("./data/raw/roads/GRIP4_GlobalRoads.gdb")
roads <- roads[roads$GP_RTP == 1, ]


saveRDS(dist_roads, "./data/output/distance_roads.RDS")

saveRDS(roads, "./data/output/motorways_sf.RDS")
raster_roads <- rasterize(roads[0],
                          chelsa_10km,
                          fun = 'count')
seq_test <- round(seq(1, 573396, length = 11))
for(i in 1:10)
{
  cat(paste0(Sys.time(), " -- iteration ", i, "\n"))
  raster_roads_test <- rasterize(roads[seq_test[i:(i + 1)]][0],
                                 chelsa_10km,
                                 fun = function(x, ...) 1)
  png(paste0("./data/raw/roads/roads_rasterized_part", i, ".png", h = 1080, w = 1920))
  plot(raster_roads_test)
  dev.off()
}


chelsa_100km <- aggregate(chelsa_10km, fact = 10)

a <- roads[1:100, ]

raster_roads_test <- rasterize(roads[1],
                               chelsa_10km,
                               fun = sum)

raster_roads2 <- raster_roads_test
raster_roads2[raster_roads2 > 1] <- 1
plot(raster_roads2)

dist_roads <- distance(raster_roads2)
saveRDS(dist_airports,
        "./data/output/distance_to_airports.RDS")


