rm(list = ls())

#library(raster)
#devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
library(biomod2)
library(tidyverse)
library(ggplot2)
library(Rarity)
library(terra)

# file path
sdm_path <- "Z:/R/sdmrisk/"

##### Data preparation ######

##### open data

baseline = rast(paste0(sdm_path, "data/output/baseline_scale_moll.tif"))
var <- readRDS(paste0(sdm_path, "data/output/var_names.RDS"))
names(baseline) <- var
sp_info <- readRDS(paste0(sdm_path, "data/sp_list_to_model.RDS"))
bg_points <- readRDS(paste0(sdm_path, "data/bg_points_bias.RDS"))
groups <- readRDS(paste0(sdm_path, "data/collinearity_groups.RDS"))


##### Select models

#models = c("GLM", "GBM", "GAM", "MARS", "MAXENT")
models = c("GLM", "GBM", "GAM", "MAXENT")

##### calibrate background points 

bg_cells <-   cellFromXY(baseline,
                         bg_points)
bg_data <- as.data.frame(bg_points)
bg_data$Observed <- 0
plot(bg_data$x, bg_data$y)
dev.off()

##### Select one var per intercorrelated group

# on garde toutes les variables seules
# on sélectionne (arbitrairement) une var par groupe de var intercorrélées
# voir cahier (29/03/23) pour le choix des var => on prend celle qui traduit un max/min 
# plutôt que de prendre des var moyennes

# On met dans un objet tous les groupes qui ont plusieurs variables
sel_groups <- groups[which(sapply(groups, length) > 1)]
sel_groups
# On prépare un tableau qui va contenir les variables sélectionnées pour chaque groupe
# on ne fait la sélection que sur les espèces avec assez de données
# voir ce qu'on fait ensuite pour les autres ?
kept_vars <- data.frame(species = sp_info$sp_key[which(sp_info$occ_qty != "Insufficient")])

kept_vars$group1 <- "bio6"
kept_vars$group2 <- "bio4"
kept_vars$group3 <- "bio5"
kept_vars$group4 <- "bio16"
kept_vars$group5 <- "bio15"
kept_vars$group5 <- "pop"

# On récupère maintenant toutes les variables pas corrélées aux autres
lone_vars <- unlist(groups[which(sapply(groups, length) == 1)])

# Et on va les rajouter aux variables sélectionnées dans chaque groupe de variables intercorrélées
if(nrow(kept_vars) == 1){
  kept_vars[, lone_vars] <- lone_vars
  vars_per_sp <- kept_vars
  } else {
  vars_per_sp <- data.frame(kept_vars,
                            sapply(lone_vars, rep,  nrow(kept_vars)))
  }

# remove target group from variable selection
vars_per_sp$plant <- vars_per_sp$terr_vert <- vars_per_sp$invert_fish <- NULL

saveRDS(vars_per_sp, file = paste0(sdm_path, "data/vars_per_sp.RDS"))



##### Calibration of one model per sp with all selected var #####

vars_per_sp <- readRDS(paste0(sdm_path, "data/vars_per_sp.RDS"))


for (sp in vars_per_sp$species[41:50]){
  # sp = "1311649"
  sp_ok = paste0("sp", sp)
  
  # Load sp occurrences
  sp_occ <- readRDS(paste0(sdm_path, "data/occ_rast/occ_spk_", sp))  
  # keep all not intercorrelated variables
  sp_env_stack <- baseline[[
    vars_per_sp[vars_per_sp$species == sp, 
                2:ncol(vars_per_sp)]]]

  # 1. Biomod2 doit-il générer des pseudoabsences ?
  # Non on a des points de bg
  
  runs_PA <- 0
  nb_PA <- 0
  runs_CV <- 2
  
  # 2. retirer les pts de bg qui sont dans des cellules avec présence
  
  # cellules avec présence
  species_cells <- cellFromXY(baseline,
                              sp_occ[, c("x", "y")])
  
  # if(any(bg_cells %in% species_cells)){
  #   bg_safe <- bg_data[-which(bg_cells %in% species_cells), ]
  # } else {
  #   bg_safe <- bg_data
  # }
  # 
  # # nombre de cell pour chaque condition
  # prNum <- nrow(sp_occ) # number of presences
  # bgNum <- nrow(bg_safe) # number of backgrounds

  # final points to consider: bind presence and background
  F_points <- bind_rows(sp_occ, bg_data)
  
  # 3. adequate occurrence format for biomod2
  
  coorxy <- F_points[, c("x", "y")]
  P_points <- F_points[, "Observed"]
  
  
  if(!dir.exists(paste0(sdm_path, "var_selection/all_vars/", sp_ok))){
    dir.create(paste0(sdm_path, "var_selection/all_vars/", sp_ok), recursive = T)
  }
  # setwd(paste0(initial_wd, "/var_selection/all_vars/", sp))
  # 
  
  MyExpl <- sp_env_stack[[
    sample(names(sp_env_stack), 
           nlyr(sp_env_stack))]]

  run_data <- BIOMOD_FormatingData(resp.name = sp_ok,
                                   resp.var = P_points, 
                                   expl.var = MyExpl, 
                                   dir.name = paste0(sdm_path, "var_selection/all_vars"),
                                   resp.xy = coorxy,
                                   PA.nb.rep = runs_PA,
                                   PA.nb.absences = nb_PA,
                                   PA.strategy = 'random')
  
  saveRDS(run_data, file = paste0(sdm_path, "var_selection/all_vars/", sp_ok, "/run_data.RDS"))
  
  myBiomodOptions <- BIOMOD_ModelingOptions()
  
  model_runs <- BIOMOD_Modeling(bm.format = run_data, 
                                modeling.id = "1", 
                                models =  models,
                                bm.options = myBiomodOptions, # Options de modélisation à mettre dans cet argument
                                nb.rep = runs_CV, 
                                data.split.perc = 100, 
                                do.full.models = FALSE, 
                                weights = NULL, 
                                prevalence = 0.5, 
                                var.import = 10, 
                                nb.cpu = 4, 
                                do.progress = TRUE)
  
  saveRDS(model_runs, file = paste0(sdm_path, "var_selection/all_vars/", sp_ok, "/model_runs.RDS"))
  print(Sys.time())
}





##### Etape 3 : sélection des variables à importance significative #####
sel_vars<-list()

for (i in 1:length(sp_list$sp))
{
  sp <- sp_list$sp[i]
  
  sp = "9613389"
  sp_ok = paste0("sp", sp)
  
  model_runs <- readRDS(paste0(sdm_path, "var_selection/all_vars/", sp_ok, "/model_runs.RDS"))
  
  gg_varimp <- get_variables_importance(model_runs)
  
  
  colnames(gg_varimp) <- c("id", 
                           "PA.Run",
                           "CV.Run", "Model", 
                           "Variable", 
                           "VI.run",
                           "Variable.importance")
  
  gg_varimp$Variable <- reorder(gg_varimp$Variable,
                                gg_varimp$Variable.importance,
                                median,
                                na.rm=TRUE)
  
  p <-  ggplot(gg_varimp, aes(y = Variable, x = Variable.importance)) +
    geom_boxplot() + geom_jitter(alpha = .2, aes(col = Model)) +
    theme_bw() + ggtitle(sp) + scale_color_brewer(palette = "Set2")
  
  
  
  png(paste0(sdm_path, "graphs/variable_importance_", sp_ok, ".png"))
  print(p)
  dev.off()
  
  gg_varimp$Variable <- factor(gg_varimp$Variable, levels = rev(levels(gg_varimp$Variable)))
  
  # On ne garde que les variables importantes à plus de x% pour au moins 50% des modèles
  # Alternative : garder les X meilleures variables
  quantile.50 <- aggregate(Variable.importance ~ Variable, data = gg_varimp, FUN = median)
  
  # Evaluation des interactions entre variables
  var_interactions <- spread(gg_varimp, Variable, Variable.importance)
  
  png(paste0(sdm_path, "graphs/variable_interactions_", sp_ok, ".png"),
      width = 900, height = 900)
  corPlot(var_interactions[, -c(1:5)], method = "pearson")
  dev.off()  
  if(any(cor(var_interactions[, -c(1:5)], use = "na.or.complete") <= -0.4))
  {
    warning(paste0("Negative interactions among variables have been detected for
                   species ", sp, ".\nYou should check the correlation plots among
                   variable importances."))
    # Retirer les variables de manière itérative: retirer la variable avec 
    # interaction négative à l'importance la plus faible, puis recalculer les 
    # corrélations, et ainsi de suite.
  }
  
  sel_vars[[sp]] <- as.character(quantile.50$Variable[which(quantile.50$Variable.importance >= 0.1)])
}



saveRDS(sel_vars, file = paste0("data/selected_variables.RDS"))


##########################

n_permut <- 6
nb_model_repetitions <- 5

var_imp <- array(dim = c(length(sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")]),
                         n_permut * nb_model_repetitions, 
                         nlyr(baseline)),
                 dimnames = list(sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")],
                                 1:(n_permut * nb_model_repetitions),
                                 names(baseline)))



# Setup parallel computing
n.cores <- parallel::detectCores() - 2
startcluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = startcluster)

# Prepare list to save selected variables in output
selected_variables <- list()

sp_info$nb_variables <- NA
sp_info$variables <- NA
sp_info$check_var_cor <- FALSE

# Need to close all graphic devices for the loop to work
for(i in 1:10)
{
  try(dev.off())
}
importance_cutoff <- 0.05
k <- 0

pdf(paste0(sdm_path, "graphs/variable_importance.pdf"))
for(sp in sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")])
{
  cat(paste0(Sys.time(), " --- Taxon ", sp, "\n"))
  
  k <- k + 1
  message(round(k / length(sp_info$sp_key[which(sp_info$occ_qty == "Sufficient")]), 4) * 100, "%")
  
  sp_target_group <- sp_info$target_group[which(sp_info$sp_key == sp)]

  
  # Load taxon occurrences
  
  sp = "2498305"
  
  # Load taxon occurrences
  sp_occ <- readRDS(paste0(sdm_path, "data/occ_rast/occ_spk_", sp))
  
  sp_env <- terra::extract(baseline,
                       sp_occ[, c("x", "y")])
  
  species_cells <- cellFromXY(baseline,
                              sp_occ[, c("x", "y")])
  
  # retirer les pts de bg qui sont dans des cellules avec présence
  if(any(bg_cells %in% species_cells)){
    bg_env <- bg_env_all[-which(bg_cells %in% species_cells), ]
  } else {
    bg_env <- bg_env_all
  }
  
  # species_occ <- c(rep(1, nrow(taxon_env)),
  #                  rep(0, nrow(bg_env)))
  
  model_table <- rbind(data.frame(occ = 1, sp_env[, variables]),
                       data.frame(occ = 0, bg_env[, variables]))
  model_table$occ <- as.factor(model_table$occ)
  
  prNum <- as.numeric(table(model_table$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(model_table$occ)["0"]) # number of backgrounds
  
  # cwt <- c("1" = 1, "0" = prNum / bgNum)
  samsize <- c("0" = prNum, "1" = prNum)
  message("Variable importance permutations...\n")
  
  
  
  
  
  
  for(md in 1:nb_model_repetitions) {
    message(paste0("\nRun ", md, "\n"))
    mod_rf <- NULL
    # convert the response to factor for RF model to return probabilities
    
    
    if(inherits(try(
      mod_rf <- randomForest(occ ~ ., 
                             data = model_table,
                             ntree = 1000,
                             sampsize = samsize,
                             replace = TRUE)
    ), "try-error")){
      print(paste("Error for taxon", taxon))
    }
    
    base_pred <- predict(mod_rf,
                         model_table[, variables[[taxo_order]]],
                         type = "prob")
    
    
    for(var in variables[[taxo_order]])
    {
      cat(".")
      varimp <- foreach(i = 1:n_permut,
                        .combine = 'c',
                        .packages = 'randomForest') %dopar% 
        {
          model_permut <- model_table[sample(1:nrow(model_table), 10000, replace = FALSE), 
                                      variables[[taxo_order]]]
          model_permut[, var] <- sample(model_permut[, var], 
                                        replace = TRUE)
          permut_pred <- predict(mod_rf,
                                 model_permut,
                                 type = "prob")
          
          1 - cor(base_pred[rownames(model_permut), 2], permut_pred[, 2])
        }
      var_imp[taxon, (1:n_permut) + (md - 1) * n_permut, var] <- varimp
    }
  }
  message("\nVariable importance done")
  
  ggvi <- reshape2::melt(var_imp[taxon, ,])
  ggvi <- droplevels(ggvi[-which(is.na(ggvi$value)), ])
  colnames(ggvi) <- c("Permutation_run", "Variable", "Variable_importance")
  
  ggvi$Variable <- reorder(ggvi$Variable,
                           ggvi$Variable_importance,
                           median,
                           na.rm=TRUE)
  
  
  print(ggplot(ggvi) +
          geom_boxplot(aes(y = Variable, x = Variable_importance)) +
          theme_minimal() +
          ggtitle(taxon) +
          xlab("Variable importance") +
          ylab("Variable") +
          geom_vline(xintercept = importance_cutoff, linetype = 3))
  
  quantile.50 <- aggregate(Variable_importance ~ Variable, data = ggvi, FUN = median)
  quantile.50 <- quantile.50[order(quantile.50$Variable_importance, decreasing = TRUE), ]
  
  # Select variable with median importance above cutoff
  selected_variables[[taxon]] <- as.character(quantile.50$Variable[which(quantile.50$Variable_importance >= importance_cutoff)])
  
  # If there are more than 1 variable per 10 occurrence records, we reduce the nb of vars
  if(length(selected_variables[[taxon]]) > floor(taxa_list$n[which(taxa_list$taxon == taxon)] / 10))
  {
    message("Too many variables for the number of occurrences, reducing the number of variables...")
    selected_variables[[taxon]] <- 
      as.character(quantile.50$Variable[
        which(quantile.50$Variable_importance >= importance_cutoff)])[
          1:floor(taxa_list$n[which(taxa_list$taxon == taxon)] / 10)]
  }
  
  
  # ----- Variable interactions -----
  
  ggvi <- ggvi[which(ggvi$Variable %in% selected_variables[[taxon]]), ]
  
  if(length(selected_variables[[taxon]]) > 1)
  {
    var.interactions <- tidyr::spread(ggvi, Variable, Variable_importance)
    
    cors <- cor(var.interactions[, -1], use = "na.or.complete")
    cors[upper.tri(cors)] <- NA
    
    if(any(cors <= -0.4, na.rm = TRUE))
    {
      message("Negative correlations detected among selected variable: removing less important & correlated variables")
      
      png(paste0("./TDMs/graphs/variable_interactions_", taxon, ".png"),
          h = 1000, w = 1000)
      corPlot(var.interactions[, -1], method = "pearson")
      dev.off()  
      
      
      taxa_list$check_var_cor[which(taxa_list$taxon == taxon)] <- TRUE
      
      neg_cors <- which(cors <= -0.4, arr.ind = T)
      continue <- TRUE
      remove_vars <- NULL
      while(continue){
        # for(pb_cors in 1:nrow(neg_cors))
        # {
        first_cor <- which(rowSums(neg_cors) == max(rowSums(neg_cors)))
        if(length(first_cor) > 1)
        {
          red_neg_cors <- neg_cors[first_cor, , drop = FALSE]
          first_cor <- first_cor[which(red_neg_cors[, 1] == max(red_neg_cors))]
        }
        pb_vars <- c(rownames(cors)[neg_cors[first_cor, 1]],
                     colnames(cors)[neg_cors[first_cor, 2]])
        pb_vars_imp <- quantile.50$Variable_importance[which(quantile.50$Variable %in% pb_vars)]
        
        
        remove_vars <- c(remove_vars,
                         pb_vars[which(pb_vars_imp == min(pb_vars_imp))])
        # We need to re-assess if there are still negative correlations from there
        cors <- cors[-which(rownames(cors) %in% remove_vars), -which(rownames(cors) %in% remove_vars)]
        if(any(cors <= -0.4, na.rm = TRUE))
        {
          neg_cors <- which(cors <= -0.4, arr.ind = T)
        } else
        {
          continue <- FALSE
        }
      }
      remove_vars <- unique(remove_vars)
      selected_variables[[taxon]] <-  selected_variables[[taxon]][-which(selected_variables[[taxon]] %in% remove_vars)]
      
    }
  }
  
  taxa_list$nb_variables[which(taxa_list$taxon == taxon)] <- length(selected_variables[[taxon]])
  taxa_list$variables[which(taxa_list$taxon == taxon)] <- paste0(selected_variables[[taxon]], collapse = ",")
  
}

dev.off()

parallel::stopCluster(startcluster)

saveRDS(taxa_list, "./TDMs/data/taxa_list_with_vars.RDS")
saveRDS(selected_variables, "./TDMs/data/selected_variables.RDS")
saveRDS(var_imp, "./TDMs/data/variable_importance.RDS")