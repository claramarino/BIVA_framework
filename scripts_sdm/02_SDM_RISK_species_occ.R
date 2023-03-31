rm(list = ls())

library(terra)
library(raster)
library(tidyverse)
library(sf)

# file path
sdm_path <- "Z:/R/sdmrisk/"

##### open baseline

baseline = rast(paste0(sdm_path, "data/output/baseline_scale_moll.tif"))
var <- readRDS(paste0(sdm_path, "data/output/var_names.RDS"))
names(baseline) <- var
crs(baseline)

##### open all occurrences 

# species key
fold = "Output/Exposure/Exposure_raw/"
df_all <- readRDS(paste0(fold, "RISK_21_grid_cells_310_IAS_55km"))
sp_302 <- unique(df_all$new_key)

# file list occurrence points
occ_path <- "Output/Occurrences_clean_taxo_ok/"
occ_files <- list.files(occ_path)

# target group
tg <- readRDS(paste0(sdm_path, "data/output/tg/target_group_info"))
# store information related to occ number and final sp to model
df <- data.frame()

for (sp in sp_302){
  # sp = "2498305" # sp for testing the loop
  
  tg_sp <- tg$group[tg$new_key==sp]
  occ_sp <- occ_files[grepl(paste0("spk_", sp), occ_files)]
  occ <- readRDS(paste0(occ_path, occ_sp))
  
  # convert occ data frame into sf object
  occ_sf <- occ %>%
    mutate_at(vars(LONG, LAT), as.numeric) %>%   # coordinates must be numeric
    st_as_sf(
      coords = c("LONG", "LAT"),
      agr = "constant",
      crs = "EPSG:4326",
      stringsAsFactors = FALSE,
      remove = TRUE) %>% 
    # add a field for observations
    mutate(Observed = 1) %>%
    # filter for coordinate uncertainty
    filter(coordinateUncertaintyInMeters < 22500)
  
  # ggplot(occ_sf)+geom_sf()
  
  occ_sf_proj <- project(vect(occ_sf), y = crs(baseline))
  
  sp_env <- rasterize(x = occ_sf_proj, 
                      y = baseline,
                      field = "Observed",
                      fun = sum) # Bien vérifier le résultat et ajuster la fonction
  
  # Graphe illustrant le nombre de présences par pixel
  # plot(sp_env, col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"))
  
  # On ne garde qu'une seule valeur par pixel
  sp_env[sp_env > 1] <- 1
  # On ajoute le nom de l'espèce au raster
  names(sp_env) <- sp
  
  # On ajoute les données environnementales aux présences rasterisées
  sp_env <- c(sp_env, baseline)
  
  # sp_env
  
  # 2. Suppression des présences pour lesquelles on n'a pas de données environnementales
  # Récupération des coordonnées de toutes les cellules
  coorXY <- xyFromCell(baseline, 1:ncell(baseline))
  # Transformation du raster en data.frame pour tester les NAs
  sp_env_df <- values(sp_env)
  
  if(any(is.na(sp_env_df[, "bio1"]) & !is.na(sp_env_df[, sp]))){
    cat("Some points are in pixels without environmental values\n")
  }
  
  # On supprime les cellules pour lesquelles on n'a pas de données environnementales
  coorXY <- coorXY[-which(is.na(sp_env_df[, "bio1"])), ]
  sp_env_df <- sp_env_df[-which(is.na(sp_env_df[, "bio1"])), ]
  
  # Nombre de cellules résultantes :
  cat(sp, "\nNumber of pixels of presence:",
      "\n - Initial: ", nrow(occ_sf),
      "\n - After rasterisation: ", length(which(sp_env_df[, 1] == 1)), "\n")
  
  # 3. Récupération des occurrences rasterisées et écriture sur le disque
  P_points <- data.frame(
    # D'abord on récupère les coordonnées XY qui correspondent à nos cellules de présences/absences
    coorXY[which(!is.na(sp_env_df[, sp])), ],
    # Ensuite, on récupère la colonne qui indique présence/absence pour chaque cellule
    Observed = sp_env_df[which(!is.na(sp_env_df[, sp])), sp]) # On récupère les occurrences ici
  
  # on ajoute le target group pour la suite 
  if (nrow(P_points)>0){P_points %>% mutate(target_group = tg_sp)} 

  
  saveRDS(P_points, file = paste0(sdm_path, "data/occ_rast/occ_spk_", sp))
  
  # get info on species occurrences
  info <- data.frame(
    sp_key = sp,
    nb_occ = nrow(occ_sf_proj),
    nb_cells = length(which(sp_env_df[, 1] == 1)),
    target_group = tg_sp
  )
  df <- bind_rows(df, info)
  
  saveRDS(df, file = paste0(sdm_path, "data/occ_rast/info_occ_all_sp"))
  
  }


# nb of species with enough occurrences at this resolution

nrow(df %>% filter(nb_cells>=30))
nrow(df %>% filter(nb_cells>=10))


df$occ_qty <- ifelse(df$nb_cells >= 30,
                            "Sufficient",
                            ifelse(df$nb_cells >= 10, "Low",
                                          "Insufficient"))
df$occ_qty <- factor(df$occ_qty,
                            levels = c("Sufficient", "Low", "Insufficient"))
table(df$occ_qty,
      df$target_group)

saveRDS(df, file = paste0(sdm_path, "data/occ_rast/info_occ_all_sp"))
