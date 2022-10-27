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
deg = "1" # can be "1" or "01"
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
  mutate(expo_tot = SR_tot_ias_area + range_med + med_iasa_tot) %>%
  mutate(expo_tot_norm = (expo_tot/max(expo_tot))) %>%
  dplyr::select(x, y, expo_tot_norm, range_med, med_iasa_tot, SR_tot_ias_area)


# ajouter les résidues de sensibilité par rapport à la richesse tot
mod_ias_tot <- lm(SR_iasa_bmr_by_area~SR_tot_bmr_by_area, data = sensi)
mod_iast_tot <- lm(SR_iast_bmr_by_area~SR_tot_bmr_by_area, data = sensi)

plot(mod_ias_tot)
summary(mod_ias_tot)
res_sensi <- residuals.lm(mod_ias_tot)
res_sensi2 <- residuals(mod_ias_tot2)
res_sensi_t <- residuals.lm(mod_iast_tot)
str(res_sensi)
cor.test(res_sensi, res_sensi2)

# normalité des résidus ?
hist(res_sensi)

sensi_fin <- sensi %>% 
  dplyr::select(x, y, SR_tot_bmr_by_area, SR_iasa_bmr_by_area,
                SR_iast_bmr_by_area) %>%
  mutate(res_iasa_tot = res_sensi,
         res_iast_tot = res_sensi_t)

# separate by main /isl

expo_isl <- inner_join(expo_fin, isl, by=c("x","y"))
expo_main <- inner_join(expo_fin, main, by=c("x","y"))
sensi_isl <- inner_join(sensi_fin, isl, by=c("x","y"))
sensi_main <- inner_join(sensi_fin, main, by=c("x","y"))


pmain <- ggplot(data = sensi_main) +
  geom_tile(aes(x = x, y = y, fill = res_iasa_tot)) +
  scale_fill_viridis_c()+
  theme_map()

pisl <- ggplot(data = sensi_isl) +
  geom_tile(aes(x = x, y = y, fill = res_iasa_tot)) +
  scale_fill_viridis_c()+
  theme_classic()

hist(sensi_main %>% filter(SR_iast_bmr_by_area>0) %>% pull(res_iast_tot), n=50)
hist(sensi_isl %>% filter(SR_iast_bmr_by_area>0) %>% pull(res_iast_tot), n=50)
hist(sensi_main %>% pull(res_iast_tot), n=50)
hist(sensi_isl %>% pull(res_iast_tot), n=50)


sensi_all <- bind_rows(sensi_isl, sensi_main)

ggplot(sensi_all, aes(x=type, y = res_iasa_tot))+
  geom_violin() +
  geom_boxplot(alpha = .5)
t.test(res_iasa_tot~type, data = sensi_all)

# remove all cells with sensitivity = 0
ggplot(sensi_all %>% filter(SR_iasa_bmr_by_area>0), 
       aes(x=type, y = res_iasa_tot))+
  geom_violin() +
  geom_boxplot(alpha = .5)
t.test(res_iasa_tot~type, data = sensi_all %>% filter(SR_iasa_bmr_by_area>0))

ggplot(sensi_all %>% filter(SR_iast_bmr_by_area>0), 
       aes(x=type, y = res_iast_tot))+
  geom_violin() +
  geom_boxplot(alpha = .5)
t.test(res_iast_tot~type, data = sensi_all%>% filter(SR_iast_bmr_by_area>0))




rm(main, isl)
expo_all <- bind_rows(expo_isl, expo_main)


ggplot(expo_all, aes(x=type, y = expo_tot_norm))+
  geom_violin() +
  geom_boxplot(alpha = .5)
ggplot(expo_all, aes(x=type, y = range_med))+
  geom_violin() +
  geom_boxplot(alpha = .5)


t.test(expo_tot_norm~type, data = expo_all)





####### Explore differences between taxa #######



rm(main)

ggplot()+ 
  geom_smooth(data = sensi_isl %>% filter(y > -66.56333 & y < 66.56333 ), 
              aes(x = y, y = res_iasa_tot, color = "island"))+ 
  geom_smooth(data = sensi_main %>% filter(y > -66.56333 & y < 66.56333 ),
              aes(x = y, y = res_iasa_tot, color = "mainland"))+
  scale_colour_manual(name="Type",
                      values=c(island="turquoise", mainland="firebrick"))+
  coord_flip()+
  theme_classic()



######### Brouillon account for number of species #########

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
