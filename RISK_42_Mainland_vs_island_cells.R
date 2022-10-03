# separate island vs mainland

rm(list=ls())


library(tidyverse)

############## Get island vs mainland cells ##############
library(raster)
library(sf)

# open continent only
contin <- st_read("Z:/THESE/5_Data/Distribution_spatiale/Shpfiles_iles_continent/Continents_only_weigelt.shp")
plot(contin)

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

# select cells from continent

extract_cells_main <- function(poly, rast){
  # from a raster containing 0 for each cell
  # keep cells that are in mainland
  masked_rast <- mask(rast, poly)
  
  # for saving memory, save a df containing all cells with values
  df_main <- as.data.frame(masked_rast, xy = TRUE) %>%
    # remove empty cells
    dplyr::filter(!is.na(layer)) %>%
    dplyr::mutate(type = "mainland") %>%
    dplyr::select(-layer)
  
  return(df_main)
}
extract_cells_isl <- function(poly, rast){
  # from a raster containing 0 for each cell
  # keep cells that are not in mainland (isl + ocean)
  masked_rast <- mask(rast, poly)
  df_isl <- as.data.frame(masked_rast, xy = TRUE) %>%
    # keep only empty cells
    dplyr::filter(is.na(layer)) %>%
    dplyr::mutate(type ="island") %>%
    dplyr::select(-layer)
  
  return(df_isl)
}


main_cells <- extract_cells_main(contin, raster0.1)
isl_cells <- extract_cells_isl(contin, raster0.1)


main_cells1 <- extract_cells_main(contin, raster1)
isl_cells1 <- extract_cells_isl(contin, raster1)

saveRDS(main_cells, "Output/Mainland_island/RISK_42_mainland_cells_r01")
saveRDS(isl_cells, "Output/Mainland_island/RISK_42_isl_ocean_cells_r01")
saveRDS(main_cells1, "Output/Mainland_island/RISK_42_mainland_cells_r1")
saveRDS(isl_cells1, "Output/Mainland_island/RISK_42_isl_ocean_cells_r1")


############### Separate Expo/Sensi by mainland vs islands ###################

# select degree resolution
deg = "01" # can be "1" or "01"
# select the type of normalization
norm = "log" # can be "rank", "log", or "lin"

# load exposure 
expo <- readRDS(paste0("Output/Exposure/RISK_24_expo_norm_", norm,"_r", deg))
# load sensitivity
sensi <- readRDS(paste0("Output/Sensitivity/RISK_33_sensitivity_norm_", norm,"_r", deg))

# load main/isl distinction
main <- readRDS(paste0("Output/Mainland_island/RISK_42_mainland_cells_r", deg))
isl <- readRDS(paste0("Output/Mainland_island/RISK_42_isl_ocean_cells_r", deg))


expo_fin <- expo %>% ungroup() %>%
  mutate(expo_tot = SR_tot_ias + range_med + med_iasa_tot) %>%
  mutate(expo_tot_norm = (expo_tot/max(expo_tot))) %>%
  dplyr::select(x, y, expo_tot_norm)


sensi_fin <- sensi %>% dplyr::select(x, y, SR_tot_bmr, SR_iasa_bmr)

# separate by min /isl

expo_isl <- inner_join(expo_fin, isl, by=c("x","y"))
expo_main <- inner_join(expo_fin, main, by=c("x","y"))
sensi_isl <- inner_join(sensi_fin, isl, by=c("x","y"))
sensi_main <- inner_join(sensi_fin, main, by=c("x","y"))


ggplot(data = sensi_main) +
  geom_raster(aes(x = x, y = y, fill = SR_iasa_bmr)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()

ggplot(data = sensi_isl) +
  geom_raster(aes(x = x, y = y, fill = SR_iasa_bmr)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  theme_classic()


hist(sensi_fin %>% pull(SR_iasa_bmr), n=50)
head(sensi_isl %>% filter(SR_tot_bmr>0 & SR_iasa_bmr>0) %>% arrange(SR_iasa_bmr))
#%>% filter(SR_tot_bmr>0 & SR_iasa_bmr>0)
hist(sensi_main %>% filter(SR_tot_bmr>0 & SR_iasa_bmr>0) %>% pull(SR_iasa_bmr), n=50)
hist(sensi_isl %>% filter(SR_tot_bmr>0 & SR_iasa_bmr>0) %>% pull(SR_iasa_bmr), n=50)

summary(sensi_main)

par(mfrow=c(1,3))
hist(sensi_fin %>% filter(SR_tot_bmr>0) %>% pull(SR_iasa_bmr), n=30)
hist(sensi_isl %>% filter(SR_tot_bmr>0) %>% pull(SR_iasa_bmr), n=30)
hist(sensi_main %>% filter(SR_tot_bmr>0) %>% pull(SR_iasa_bmr), n=30)

hist(sensi_fin %>% filter(SR_tot_bmr>0) %>% pull(SR_tot_bmr), n=30)
hist(sensi_isl %>% filter(SR_tot_bmr>0) %>% pull(SR_tot_bmr), n=30)
hist(sensi_main %>% filter(SR_tot_bmr>0) %>% pull(SR_tot_bmr), n=30)

ggplot(sensi_main)+
  geom_point(aes(x= SR_tot_bmr, y =SR_iasa_bmr))


# account for total nb of species ?

sensi_fin_2 <- sensi_fin %>%
  mutate(SR_iasa_bmr_norm = if_else(SR_tot_bmr>0, SR_iasa_bmr/SR_tot_bmr, SR_iasa_bmr))
max_norm = max(sensi_fin_2$SR_iasa_bmr_norm)
sensi_fin_2 <- sensi_fin_2 %>%
  mutate(SR_iasa_bmr_norm2 = SR_iasa_bmr_norm/max_norm)


hist(sensi_fin_2 %>% filter(SR_tot_bmr>0 & SR_iasa_bmr>0) %>% pull(SR_iasa_bmr_norm2), n=50)
