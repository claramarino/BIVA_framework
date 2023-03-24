library(raster)

#***************************
#* -------Baseline-------- *
#***************************

# Climate variables
chelsa.names <- paste0("bio", sprintf("%02d", 1:19))
climate.vars <- paste0("bio", 1:19)

for (j in chelsa.names)
{
  if(j == chelsa.names[1]) # Pour le géoréférencement
  {
    env.stack <- stack(raster(paste0("./data/raw/bioclim/CHELSA_baseline_", j, ".tif")))
  } else
  {
    env.stack <- addLayer(env.stack, 
                          raster(paste0("./data/raw/bioclim/CHELSA_baseline_", j, ".tif")))
  }
}

NAvalue(env.stack) <- -32768
env.stack <- setMinMax(env.stack)


names(env.stack) <- climate.vars
env.stack <- aggregate(env.stack,
                       fact = 10)

writeRaster(env.stack,
            "./data/output/baseline_10km_bioclim")
env.stack <- stack("./data/output/baseline_10km_bioclim")
env.stack <- readAll(env.stack)
saveRDS(env.stack,
        "./data/output/baseline_10km_bioclim.RDS")

env.stack <- readRDS("./data/output/baseline_10km_bioclim.RDS")

# Land use variables

list.files(paste0(sdm_path, "data/output/lu/"))
lu.stack <- terra::rast()

rlist <- list()
for (i in 1:length(chunks)){
  rlist <- c(rlist, list(readRDS(chunks[i]) ))
}

lu.stack <- readRDS("./data/output/baseline_2015_lu.RDS")

# Distance to ports and airports
dist_ports <- readRDS("./data/output/distance_to_ports.RDS")
dist_airports <- readRDS("./data/output/distance_to_airports.RDS")

# Distance to roads
dist_roads <- readRDS("./data/output/distance_roads.RDS")

# Population
pop_baseline <- raster("./data/raw/population/BaseYear_1km/baseYr_total_2000.tif")
pop_baseline <- aggregate(pop_baseline,
                          fun = sum, # Unit is number of person so when we aggregate we have to SUM, not AVERAGE
                          fact = 10) 

pop_baseline <- extend(pop_baseline,
                       env.stack)
pop_baseline <- resample(pop_baseline,
                         env.stack,
                         method = "ngb")


env.stack <- stack(env.stack,
                   lu.stack,
                   lu.stack[["PFT1"]] + lu.stack[["PFT7"]],
                   dist_ports,
                   dist_airports,
                   pop_baseline,
                   dist_roads)
names(env.stack)[(nlayers(env.stack) - 4):nlayers(env.stack)] <- c("trees", "dist_ports", "dist_airports", "population", "dist_roads")


writeRaster(env.stack,
            "./data/output/baseline_total", overwrite = TRUE)
env.stack <- stack("./data/output/baseline_total")
env.stack <- readAll(env.stack)
saveRDS(env.stack,
        "./data/output/baseline_total.RDS")

initial.variables <- c(paste0("bio", 1:19),
                       "PFT1", # Needleleaf evergreen tree – temperate 
                       "PFT5", # Broadleaf evergreen tree – temperate
                       "PFT7", # Broadleaf deciduous tree – temperate
                       "PFT10", # Broadleaf deciduous shrub – temperate 
                       "trees",
                       "dist_ports",
                       "dist_airports",
                       "population",
                       "dist_roads")

env.stack <- env.stack[[initial.variables]]
env.stack <- virtualspecies::synchroniseNA(env.stack)
saveRDS(env.stack,
        "./data/output/baseline.RDS")

# Need to add projection



# Adding distance to roads which were calculated later on and projecting to mollweide
# Renaming land use variables
lu.vars2 <- c("net_tem", # Needleleaf evergreen tree – temperate 
              "bet_tem", # Broadleaf evergreen tree – temperate
              "bdt_tem", # Broadleaf deciduous tree – temperate
              "bds_tem") # Broadleaf deciduous shrub – temperate 

lu.vars <- c("PFT1", # Needleleaf evergreen tree – temperate 
             "PFT5", # Broadleaf evergreen tree – temperate
             "PFT7", # Broadleaf deciduous tree – temperate
             "PFT10") # Broadleaf deciduous shrub – temperate 

for (scenario in scenarios)
{
  for (gcm in gcms)
  {
    for(year in periods)
    {
      cat(paste0(Sys.time(), " - Scenario ", scenario,
                 " - ", gcm,
                 " - ", year),"\n")
      
      
      future.stack <- readRDS(paste0("./data/output/future_",
                                     scenario, "_", gcm, "_", year, ".RDS"))
      
      # future.stack <- stack(future.stack,
      #                       env.stack[["dist_roads"]])
      # names(future.stack)[nlayers(future.stack)] <- "dist_roads"
      # future.stack <- virtualspecies::synchroniseNA(future.stack)
      
      crs(future.stack) <- crs(env.stack)
      
      names(future.stack)[names(future.stack) %in% lu.vars] <- lu.vars2
      
      writeRaster(future.stack,
                  paste0("./data/output/final/future_",
                         scenario, "_", gcm, "_", year), overwrite = TRUE)
      future.stack <- readAll(future.stack)
      future.stack <- stack(future.stack)
      saveRDS(future.stack,
              paste0("./data/output/final/future_",
                     scenario, "_", gcm, "_", year, ".RDS"))
      
      future.stack.moll <- projectRaster(from = future.stack,
                                         crs = "+proj=moll",
                                         res = 10000)
      writeRaster(future.stack.moll,
                  paste0("./data/output/final/future_moll_",
                         scenario, "_", gcm, "_", year), overwrite = TRUE)
      future.stack.moll <- readAll(future.stack.moll)
      future.stack.moll <- stack(future.stack.moll)
      saveRDS(future.stack.moll,
              paste0("./data/output/final/future_moll_",
                     scenario, "_", gcm, "_", year, ".RDS"))
      
    }
  }
}


# 
# # Recreate baseline for land use
# library(data.table)
# # basecsv <- fread("./data/raw/lu/baselayerdata_region_basin_0.05deg.csv")
# basecsv <- fread("./data/raw/lu/landcover_2005_timestep.csv")
# 
# 
# lu.vars2 <- c("net_tem", # Needleleaf evergreen tree – temperate 
#               "bet_tem", # Broadleaf evergreen tree – temperate
#               "bdt_tem", # Broadleaf deciduous tree – temperate
#               "bds_tem") # Broadleaf deciduous shrub – temperate 
# 
# basecsv <- as.data.frame(basecsv)
# 
# colnames(basecsv)[1:2] <- c("y", "x")
# 
# baselu <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90, resolution = .05)
# 
# baselu2 <- rasterize(basecsv[, c("x", "y")],
#                      baselu,
#                      field = basecsv[, lu.vars2])
# baselu2 <- baselu2 * 100
# 
# 
# baselu2 <- resample(baselu2,
#                     env.stack)
# 
# initial.variables <- c(paste0("bio", 1:19),
#                        "PFT1", # Needleleaf evergreen tree – temperate 
#                        "PFT5", # Broadleaf evergreen tree – temperate
#                        "PFT7", # Broadleaf deciduous tree – temperate
#                        "PFT10", # Broadleaf deciduous shrub – temperate 
#                        "trees",
#                        "dist_ports",
#                        "dist_airports",
#                        "population",
#                        "dist_roads")
# 
# env.stack <- readRDS("./data/output/baseline.RDS")
# env.stack <- stack(env.stack[[paste0("bio", 1:19)]],
#                    baselu2,
#                    env.stack[[c("dist_ports",
#                                 "dist_airports",
#                                 "population",
#                                 "dist_roads")]])
# nmaxlayers <- nlayers(env.stack)
# future.stack <- readRDS("./data/output/final/future_rcp26_gfdl_2050.RDS")
# tmp <- future.stack[[1]]
# names(tmp) <- "tmp"
# env.stack <- addLayer(env.stack,
#                       tmp)
# env.stack <- virtualspecies::synchroniseNA(env.stack)
# env.stack <- env.stack[[1:nmaxlayers]]
# 
# writeRaster(env.stack,
#             paste0("./data/output/final/baseline"), overwrite = TRUE)
# 
# 
# env.stack <- readAll(env.stack)
# env.stack <- stack(env.stack)
# saveRDS(env.stack,
#         paste0("./data/output/final/baseline.RDS"))
# 
# env.stack.moll <- projectRaster(from = env.stack,
#                                 crs = "+proj=moll",
#                                 res = 10000)
# writeRaster(env.stack.moll,
#             paste0("./data/output/final/baseline_moll"), overwrite = TRUE)
# 
# env.stack.moll <- readAll(env.stack.moll)
# env.stack.moll <- stack(env.stack.moll)
# saveRDS(env.stack.moll,
#         paste0("./data/output/final/baseline_moll.RDS"))
# 
# 
# a <- readRDS("./data/output/final/future_moll_rcp26_gfdl_2050.RDS")
# env.stack.moll <- stack("./data/output/final/future_moll_rcp26_gfdl_2050.grd")
# env.stack.moll
# plot(env.stack.moll[[1]])
# 
# plot(future.stack.moll)