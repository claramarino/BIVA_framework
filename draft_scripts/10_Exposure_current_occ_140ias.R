# from alien occurrences and data, compute global exposure to 140 IAS

rm(list=ls())

library(tidyverse)
library(readr)
library(sp)
library(raster)
library(pbapply)
library(ggcorrplot)

#### Load alien data #####

## load alien occurrences
# for the moment, DASCO
occ_dasco <- readRDS("Output/zOld/focus_ias_range/01_ias129_occ_dasco_alien_range")


## load alien attributes
# number of native species associated
nb_nat <- readRDS("Output/focus_ias_range/01_1_nb_assoc_nat_140ias")
# total number of occurrences in gbif?
ias140_gbif <- readRDS("Output/focus_ias_range/01_2_nb_occ_140_ias")
# trophic guild?
# other attributes?

#### Bind alien attributes in one df ####

nb_occ_dasco <- occ_dasco@data %>% 
  group_by(taxon) %>% 
  summarise(nb_occ_alien = n())
alien_attrib <- left_join(nb_nat %>% rename(taxon=species),
                          left_join(ias140_gbif, nb_occ_dasco, by="taxon"),
                          by = "taxon")


#### Extract alien presences for each cell ####


# define extract cells function

# INPUT : a spatial point object and a raster of desired resolution
# OUTPUT : a df containing all the cells with a point, and with nb of pts per cell
extract_cells_pts <- function(pts, rast){
  
  # extract cells with a point inside
  c_pts <- as.data.frame(extract(rast, pts, cellnumbers=TRUE))
  # for each cell, count nb of points inside (density)
  c_pts_dens <- c_pts %>% group_by(cells) %>% summarize(n_pts = n())
  
  # set up species name
  binomial_df <- data.frame(
    binomial = rep(unique(pts$taxon), nrow(c_pts_dens)))
  
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


# initialize raster grid
raster0.1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    resolution = 0.1, vals=NULL)
raster0.1[] <- 0

# set up species list from alien occurrence data
taxon <- as.list(unique(occ_dasco$taxon))

# extract cell values for all species 
df_129_ias <- pblapply(taxon, function(x){
  pts_sp <- subset(occ_dasco, taxon == x)
  cells <- extract_cells_pts(pts_sp, raster0.1)
  return(cells)
}) # takes 45sec


# combine df for each sp into one df
df_all <- bind_rows(df_129_ias) 

str(df_all)


# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(binomial) %>%
  summarise(n_cells_range = n())

# bind with alien attributes
alien_attrib_tot <- left_join(tot_range, 
                              alien_attrib %>% rename(binomial = taxon), 
                              by = "binomial")
cor_attr <- round(cor(alien_attrib_tot %>% dplyr::select(-binomial)),4)

ggcorrplot(cor_attr)
ggplot(alien_attrib_tot, aes(x = n_cells_range, y = nb_nat)) +
  geom_point()


df_all_attrib <- left_join(df_all,
                           alien_attrib_tot,
                           by = "binomial")



# calculate cell related metrics
ias_all_agg <- df_all_attrib %>%
  group_by(x, y, cells) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_pts), # nb occ per cell (all ias sp)
    nb_occ_med = median(n_pts), # median nb occ per cell
    nb_nat_tot = sum(nb_nat),
    nb_nat_med = median(nb_nat),
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range)
    ) 

str(ias_all_agg)
summary(ias_all_agg)

# visualize correlations between all metrics
cormat <- round(cor(ias_all_agg),4)
ggcorrplot(cormat)




# plot species richness
SR_ias<- ggplot(data = ias_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$binomial))),
       x = "Longitude", y = "Latitude")
SR_ias




