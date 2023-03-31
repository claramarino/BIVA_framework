# préparaton de la block CV
rm(list = ls())

library(terra)
library(raster)
library(tidyverse)
library(sf)
library(virtualspecies)
library(blockCV)
library(dismo)


# file path
sdm_path <- "Z:/R/sdmrisk/"

##### open baseline

baseline = stack(rast(paste0(sdm_path, "data/output/baseline_scale_moll.tif")))
var <- readRDS(paste0(sdm_path, "data/output/var_names.RDS"))
names(baseline) <- var
sp_info <- readRDS(paste0(sdm_path, "data/occ_rast/info_occ_all_sp"))

plot(baseline[[1]])

##### remove collinear var

collinearity_groups <- virtualspecies::removeCollinearity(
  baseline,
  plot = TRUE,
  method= "spearman",
  multicollinearity.cutoff = 0.7)

saveRDS(collinearity_groups,paste0(sdm_path, "data/collinearity_groups.RDS"))

collinearity_groups

var_interest <- c(
  "bio6", "bio2","bio4","bio5","bio8","bio16","bio15","bio19", 
  "bare","pop","forest","grass","moss_lichen","shrubs","water",
  "invert_fish","plant","terr_vert"
)


##### Génération des blocs de cross-validation pour chaque espèce

testblocks <-  cv_spatial_autocor(
  r = rast(baseline[[var_interest]]), # raster file
  num_sample = 5000, # number of cells to be used
  plot = TRUE)


minblocksize <- min(testblocks$range_table$range)
medianblocksize <- median(testblocks$range_table$range)


bg_points <- dismo::randomPoints(baseline,
                                 n = 30000)

saveRDS(bg_points,paste0(sdm_path, "data/bg_points.RDS"))


bg_data <- sf::st_as_sf(data.frame(bg_points),
                        coords = c("x", "y"),
                        crs = raster::crs(baseline))
bg_data$Observed <- 0
plot(bg_data, pch=".")

# Choisir le nombre de folds de cross-validation
nb.of.folds = 3
unbalanced_species <- NULL
sp_info$blockCV <- "No block CV"

# Need to close all graphic devices for the loop to work
for(i in 1:10){ try(dev.off())}

pdf(paste0(sdm_path, "graphs/blocks.pdf"))
dev.new()


# Définition des folds
# uniquement pour les sp avec un nombre suffisant de cellules (n>30)

for(sp in sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")]){
  
  # sp = "2498305"
  
  # Load taxon occurrences
  sp_occ <- readRDS(paste0(sdm_path, "data/occ_rast/occ_spk_", sp))
  
  # Transformation des occurrences en objet spatial (format sf)
  pres_data <- sf::st_as_sf(sp_occ[, c("x", "y")],
                            coords = c("x", "y"),
                            crs = raster::crs(baseline))
  pres_data$Observed <- 1
  
  pa_data <- rbind(pres_data,
                   bg_data)
  
  saveRDS(pa_data,
          paste0(sdm_path, "data/blockCV/occurrences/spk_", sp, ".RDS"))
  
  dev.set(3)
  # Si pas de pseudoabsences, dans ce cas on ne crée les folds qu'une fois
  blockrun <-  try(cur.blocks <- cv_spatial(
    x = pa_data, # Données présence-absence
    column = "Observed", # Colonne contenant l'info présence/absence
    hexagon = T,
    # rasterLayer = baseline[["bio1"]], # Facultatif, quel fond de carte utiliser pour le graphe
    size = minblocksize, # Taille des blocs
    k = nb.of.folds, # Nombre de folds
    selection = "random", # Méthode d'attribution des blocs en folds. Cf. ?spatialBlock
    iteration = 100, # Nombre de runs pour répartir les blocs en folds
    biomod2 = TRUE,  # Créer une table au format biomod2 ici
    plot = TRUE))
  
  
  dev.set(2)
  # To explore later: optimize the block size
  # if(inherits(blockrun, "try-error"))
  # {
  #   blockmax <- medianblocksize
  #   blockmin <- minblocksize
  #   direction <- "reduce"
  #   ntry <- 1
  #   while(inherits(blockrun, "try-error") & ntry <= 10)
  #   {
  #     blocktest <- blockmin + (blockmax - blockmin) / 2
  #     blockrun <- try(cur.blocks <- spatialBlock(speciesData = pa_data, # Données présence-absence
  #                                                species = "Observed", # Colonne contenant l'info présence/absence
  #                                                rasterLayer = baseline[["BO21_tempmean_bdmean"]], # Facultatif, quel fond de carte utiliser pour le graphe
  #                                                theRange = blocktest, # Taille des blocs
  #                                                k = nb.of.folds, # Nombre de folds
  #                                                selection = "random", # Méthode d'attribution des blocs en folds. Cf. ?spatialBlock
  #                                                iteration = 100, # Nombre de runs pour répartir les blocs en folds
  #                                                biomod2Format = TRUE)) # Créer une table au format biomod2 ici
  #     if(inherits(blockrun, "try-error"))
  #     {
  #       blockmax <- blocktest
  #     } else
  #     {
  #       blockmin <- blocktest
  #       class(blockrun) <- append(class(blockrun), "try-error")
  #     }
  #     print(blocktest)
  #     ntry <- ntry + 1
  #   }
  # }
  
  cat("----   ", sp, "\n")
  cat("Max variation from median training presences: ")
  cat(round(max((blockrun$records$train_1 - median(blockrun$records$train_1)) / median(blockrun$records$train_1)), 2))
  
  
  # Plots disabled
  # print(blockrun$plots + ggtitle(sp))
  # str(blockrun)
  # cur.blocks.smallsize <- blockrun[-which(names(blockrun) == "plots")]

  saveRDS(blockrun, paste0(sdm_path, "data/blockCV/spk", sp, "_blocks.RDS"))
  
}
dev.off()

########################## reprendre ici

sp_info$blockCV <- "No block CV"
for(sp in sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")])
{
  blockrun <- readRDS(paste0(sdm_path, "data/blockCV/spk", sp, "_blocks.RDS"))
  
  
  cat("----   ", sp, "\n")
  print(blockrun$records)
  
  cat("Max variation from median training presences: \n")
  cat(round(max(abs(blockrun$records$train_1 - median(blockrun$records$train_1))) / median(blockrun$records$train_1), 2), "\n")
  
  if(any(abs(blockrun$records$train_1 - median(blockrun$records$train_1)) / median(blockrun$records$train_1) > 0.3))
  {
    sp_info$blockCV[which(sp_info$sp_key == sp)] <- "Unbalanced"
  } else
  {
    sp_info$blockCV[which(sp_info$sp_key == sp)] <- "Balanced"
  }
  
  cat("Average ratio train / test\n")
  cat(round(mean(blockrun$records$train_1 / rowSums(blockrun$records[, c("train_1", "test_1")])), 2), "\n")
  
  if(any(abs(blockrun$records$train_1 / rowSums(blockrun$records[, c("train_1", "test_1")]) - (1 - 1/nb.of.folds)) > 0.1))
  {
    sp_info$blockCV[which(sp_info$sp_key == sp)] <- "Unbalanced"
  } else
  {
    sp_info$blockCV[which(sp_info$sp_key == sp)] <- "Balanced"
  }
  
  
  # For those taxa where there are less than 10 evaluation points for a fold, 
  # we will make no block CV
  if(any(blockrun$records$test_1 < 10)) {
    sp_info$occ_qty[which(sp_info$sp_key == sp)] <- "Low"
  }
}

sp_info[which(sp_info$blockCV == "Unbalanced"), ]
table(sp_info$blockCV, sp_info$occ_qty)

# Print tables for unbalanced taxa
for(sp in sp_info$sp_key[which(sp_info$blockCV == "Unbalanced" &
                                 sp_info$occ_qty == "Low")])
{
  blockrun <- readRDS(paste0(sdm_path, "data/blockCV/spk", sp, "_blocks.RDS"))
  
  
  cat("----   ", sp, "\n")
  print(blockrun$records)
  
  cat("Max variation from median training presences: \n")
  cat(round(max(abs(blockrun$records$train_1 - median(blockrun$records$train_1))) / median(blockrun$records$train_1), 2), "\n")
  
  cat("Ratios train / test\n")
  cat(round(blockrun$records$train_1 / rowSums(blockrun$records[, c("train_1", "test_1")]), 2), "\n")
  
}

saveRDS(sp_info, paste0(sdm_path, "data/sp_list_to_model.RDS"))

