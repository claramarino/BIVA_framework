rm(list = ls())

library(tidyverse)
library(terra)

# file path
sdm_path <- "Z:/R/sdmrisk/"


# Climate variables

bioclim_files <- list.files(paste0(sdm_path, "data/raw/bioclim/"), full.names = T)
bio.stack <- terra::rast(bioclim_files)

names(bio.stack) <- paste0("bio", 1:19)
bio.stack
# aggregate bioclim data to 50 km
# specify na.rm=T for not removing too many cells (inclding islands...)
bio.stack50 <- terra::aggregate(bio.stack, fact = 50, na.rm=T)
bio.stack50

# Land use variables

lu_names <- c("bare", "built", "crop", "forest", "grass", 
              "moss_lichen", "shrubs", "water")
files_lu <- c()
for (lu in lu_names){
  files <- list.files(paste0(sdm_path, "data/output/lu/"), full.names = T)
  files_lu <- c(files_lu, files[grepl(lu, files)])
  }

lu.stack <- terra::rast(files_lu)
names(lu.stack) <- lu_names
lu.stack

# resample at climate resolution
lu.stack50 <- resample(lu.stack, bio.stack50[[1]], method = "bilinear")
lu.stack50


# roads

roads <- terra::rast(paste0(
  sdm_path, "data/raw/roads/grip4_total_dens_m_km2.asc"))
# resample as bioclim
roads
roads50 <- resample(roads, bio.stack50[[1]], method = "bilinear")
plot(roads)
plot(roads50, col = c("#E34A33", "#B30000"))
names(roads50) <- "roads"

# population

pop <- terra::rast(paste0(
  sdm_path, "data/raw/pop/gpw_v4_population_count_rev11_2015_2pt5_min.tif"))
pop
# aggregate as bioclim
pop50 <- aggregate(pop, 10, fun ="sum", na.rm=T) # total pop (not density) so sum needed
pop50 <- resample(pop50, bio.stack50[[1]], method = "bilinear")
names(pop50) <- "pop"
plot(pop50, col = c("#E34A33", "#B30000"))


# target group

tg_names <- c("SR_tot_plant", "SR_tot_vert_terr", "SR_fish_invert")
tg_new <- c("plant","terr_vert","invert_fish")
files_tg <- c()
for (tg in tg_names){
  files <- list.files(paste0(sdm_path, "data/output/tg/"), full.names = T)
  files_tg <- c(files_tg, files[grepl(tg, files)])
}

tg.stack <- terra::rast(files_tg)
names(tg.stack) <- tg_new

tg.stack[is.na(tg.stack)] <- 0


#### stack all environmental layers together #####
ext(bio.stack50)==ext(pop50)

env.stack <- c(bio.stack50, 
               lu.stack50, 
               # roads50, # retirer pcq bcp + de NA que les autres
               pop50, 
               tg.stack)

names(env.stack)

writeRaster(env.stack, 
            paste0(sdm_path, "data/output/baseline_envt_stack_50km.tif"),
            overwrite = T)


# synchronize NA
library(virtualspecies)
env.stack <- stack(
  rast(paste0(sdm_path, "data/output/baseline_envt_stack_50km.tif")))
plot(env.stack[["plant"]])

var <- names(env.stack)
saveRDS(var, paste0(sdm_path, "data/output/var_names.RDS"))

env.stack <- synchroniseNA(env.stack)

names(env.stack)
saveRDS(env.stack,paste0(sdm_path, "data/output/baseline.RDS"))
writeRaster(env.stack, paste0(sdm_path, "data/output/baseline.tif"),
            overwrite = T)



# Need to add projection

moll <- crs("+proj=moll")
moll

env.stack <- rast(paste0(sdm_path, "data/output/baseline.tif"))
names(env.stack) <- var

env.stack.moll <- terra::project(env.stack, moll, res = 50000)

plot(env.stack.moll[["pop"]], col = c("#E34A33", "#B30000"))

writeRaster(env.stack.moll,
            paste0(sdm_path, "data/output/baseline_moll.tif"), overwrite = TRUE)


# scale all variables from baseline
scale_params <- data.frame(mean = cellStats(stack(env.stack.moll), 'mean'),
                           sd = cellStats(stack(env.stack.moll), 'sd'))
saveRDS(scale_params, paste0(sdm_path, "data/output/scale_param_moll.RDS"))

scale_params <- data.frame(mean = cellStats(stack(env.stack), 'mean'),
                           sd = cellStats(stack(env.stack), 'sd'))
saveRDS(scale_params, paste0(sdm_path, "data/output/scale_param.RDS"))


# Scaling
baseline <- scale(stack(env.stack))
names(baseline) <- var
baseline_moll <- scale(stack(env.stack.moll))
names(baseline_moll) <- var

writeRaster(baseline,
            paste0(sdm_path, "data/output/baseline_scale.tif"), overwrite = TRUE)

writeRaster(baseline_moll,
            paste0(sdm_path, "data/output/baseline_scale_moll.tif"), overwrite = TRUE)
