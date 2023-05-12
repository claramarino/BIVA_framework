rm(list = ls())

#library(raster)
#devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
library(biomod2)
library(tidyverse)
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

# bias in bg or not ?
bg_points <- readRDS(paste0(sdm_path, "data/bg_points_bias.RDS"))
bg_points <- readRDS(paste0(sdm_path, "data/bg_points.RDS"))

groups <- readRDS(paste0(sdm_path, "data/collinearity_groups.RDS"))


##### Select models

models = c("GLM", "GBM", "GAM", "MAXENT")
# run on several computers for time saving
# models = c("GLM", "GBM", "GAM") # for virtual machines
# models = c("MAXENT") # for personnal PC (doesn't work on VM) 

##### calibrate background points 

bg_cells <-   cellFromXY(baseline,
                         bg_points)
bg_data <- as.data.frame(bg_points)
bg_data$Observed <- 0
plot(bg_data$x, bg_data$y)
#dev.off()

##### Select one var per intercorrelated group #####

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
kept_vars$group6 <- "pop"
kept_vars$terr_vert <- "terr_vert"

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

# save with target group
saveRDS(vars_per_sp, file = paste0(sdm_path, "data/vars_per_sp_tg.RDS"))

# remove target group from variable selection
vars_per_sp$plant <- vars_per_sp$terr_vert <- vars_per_sp$invert_fish <- NULL
saveRDS(vars_per_sp, file = paste0(sdm_path, "data/vars_per_sp.RDS"))



##### Calibration of one model per sp with all selected var #####

vars_per_sp <- readRDS(paste0(sdm_path, "data/vars_per_sp_tg.RDS"))


fourmi_bourdon <- c("1314773", "1340503")



for (sp in fourmi_bourdon){
  # sp = "1311649"
  sp_ok = paste0("sp", sp)
  
  # Load sp occurrences
  sp_occ <- readRDS(paste0(sdm_path, "data/occ_rast/occ_spk_", sp))  
  # keep all not intercorrelated variables
  sp_env_stack <- baseline[[
    vars_per_sp[vars_per_sp$species == sp, 
                2:ncol(vars_per_sp)]]]
  
  tg = sp_info$target_group[sp_info$sp_key==sp]
  
  gp = c("plant","terr_vert","invert_fish")
  
  to_rm <- gp[which(gp!=tg)]
  for(g in to_rm){
    sp_env_stack[[g]] <- NULL
  }

  # 1. Biomod2 doit-il générer des pseudoabsences ?
  # Non on a des points de bg
  
  runs_PA <- 0
  nb_PA <- 0
  runs_CV <- 2


  # final points to consider: bind presence and background
  F_points <- bind_rows(sp_occ, bg_data)
  
  # 3. adequate occurrence format for biomod2
  
  coorxy <- F_points[, c("x", "y")]
  P_points <- F_points[, "Observed"]
  
  
  if(!dir.exists(paste0(sdm_path, "fourmi_bourdon/", sp_ok))){
    dir.create(paste0(sdm_path, "fourmi_bourdon/", sp_ok), recursive = T)
  }
  # setwd(paste0(initial_wd, "/var_selection/all_vars/", sp))
  # 
  
  MyExpl <- sp_env_stack[[
    sample(names(sp_env_stack), 
           nlyr(sp_env_stack))]]

  run_data <- BIOMOD_FormatingData(resp.name = sp_ok,
                                   resp.var = P_points, 
                                   expl.var = MyExpl, 
                                   dir.name = paste0(sdm_path, "fourmi_bourdon"),
                                   resp.xy = coorxy,
                                   PA.nb.rep = runs_PA,
                                   PA.nb.absences = nb_PA,
                                   PA.strategy = 'random')
  
  saveRDS(run_data, file = paste0(sdm_path, "fourmi_bourdon/", sp_ok, "/run_data.RDS"))
  
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
  
  saveRDS(model_runs, file = paste0(sdm_path, "fourmi_bourdon/", sp_ok, "/model_runs.RDS"))
  print(Sys.time())
}




##### Check computation time #####

# voir le temps mis par espèce
# lié au nombre d'occurrences ?
a <- fs::dir_info(paste0(sdm_path, "fourmi_bourdon/"))
a$duree <- a$change_time - a$birth_time
for (i in 1:nrow(sp_info)){
  sp = paste0("sp",sp_info$sp_key[i])
  if(sum(grepl(sp, a$path))==0){
    sp_info$time_calc[i] <- NA
  } else {
    sp_info$time_calc[i] <- a$duree[grep(sp, a$path)]
  }
  
}

colnames(sp_info)
ggplot(sp_info)+
  geom_point(aes(x = nb_cells, y = time_calc))


hist(sp_info$nb_cells)
dev.off()



##### Etape 3 : sélection des variables à importance significative #####
a <- fs::dir_info(paste0(sdm_path, "fourmi_bourdon/"))
sp_vars <- c()
for (i in 1:nrow(sp_info)){
  sp = paste0("sp",sp_info$sp_key[i])
  if(sum(grepl(sp, a$path))>0){
    sp_vars <- c(sp_vars, sp)
  }
}

sel_vars<-list()
sp_info$nb_variables <- NA
sp_info$variables <- NA
sp_info$check_var_cor <- FALSE


# threshold for variable importance value
vimp = 0.05 # 0.1 too high?
# or nb or var to keep 
nvar = 6

for (i in 1:length(sp_vars)){
  print(i)

  sp_ok = sp_vars[i]
  sp = substring(sp_ok, 3)
  
  model_runs <- readRDS(paste0(sdm_path, "fourmi_bourdon/", sp_ok, "/model_runs.RDS"))

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
    theme_bw() + ggtitle(sp_ok) + scale_color_brewer(palette = "Set2")
  
  
  
  png(paste0(sdm_path, "fourmi_bourdon/var_imp_", sp_ok, ".png"))
  print(p)
  dev.off()
  
  gg_varimp$Variable <- factor(gg_varimp$Variable, levels = rev(levels(gg_varimp$Variable)))
  
  # On ne garde que les variables importantes à plus de x% pour au moins 50% des modèles
  # Alternative : garder les X meilleures variables
  quantile.50 <- aggregate(Variable.importance ~ Variable, data = gg_varimp, FUN = median)
  quantile.50 <- quantile.50[order(quantile.50$Variable.importance, decreasing = TRUE), ]
  
  # Select variable with median importance above cutoff (vimp)
  sel_vars[[sp_ok]] <- as.character(quantile.50$Variable[which(quantile.50$Variable.importance >= vimp)])
  
  # si trop peu de var sélectionnées, prend les 6 plus importantes
  if(length(sel_vars[[sp_ok]]) < 6){
    sel_vars[[sp_ok]] <- as.character(quantile.50$Variable[1:nvar])
  }
  
  # If there are more than 1 variable per 10 occurrence records, we reduce the nb of vars
  if(length(sel_vars[[sp_ok]]) > floor(sp_info$nb_cells[which(sp_info$sp_key == sp)] / 10))
  {
    message("Too many variables for the number of occurrences, reducing the number of variables...")
    sel_vars[[sp_ok]] <- 
      as.character(quantile.50$Variable[
        which(quantile.50$Variable.importance >= vimp)])[
          1:floor(sp_info$nb_cells[which(sp_info$sp_key == sp)] / 10)]
  }
  
  # Evaluation des interactions entre variables
  
  ggvi <- gg_varimp[which(gg_varimp$Variable %in% sel_vars[[sp_ok]]), ]
  
  if(length(sel_vars[[sp_ok]]) > 1)
  {
    var.interactions <- tidyr::spread(ggvi, Variable, Variable.importance)
    
    cors <- cor(var.interactions[, -c(1:5)], use = "na.or.complete")
    cors[upper.tri(cors)] <- NA
    
    if(any(cors <= -0.4, na.rm = TRUE))
    {
      message("Negative correlations detected among selected variable: removing less important & correlated variables")
      
      png(paste0(sdm_path, "fourmi_bourdon/variable_interaction_", sp_ok, ".png"))
      corPlot(var.interactions[, -c(1:5)], method = "pearson")
      dev.off()  
      
      sp_info$check_var_cor[which(sp_info$sp_key == sp)] <- TRUE
      
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
        pb_vars_imp <- quantile.50$Variable.importance[which(quantile.50$Variable %in% pb_vars)]
        names(pb_vars_imp) <- quantile.50$Variable[which(quantile.50$Variable %in% pb_vars)]
    
        remove_vars <- c(remove_vars,
                         names(which(pb_vars_imp == min(pb_vars_imp))))
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
      sel_vars[[sp_ok]] <-  sel_vars[[sp_ok]][-which(sel_vars[[sp_ok]] %in% remove_vars)]
      
    }
  }

  sp_info$nb_variables[which(sp_info$sp_key == sp)] <- length(sel_vars[[sp_ok]])
  sp_info$variables[which(sp_info$sp_key == sp)] <- paste0(sel_vars[[sp_ok]], collapse = ",")
}


# classic vimp <.1
saveRDS(sel_vars, file = paste0(sdm_path, "fourmi_bourdon/selected_variables.RDS"))
saveRDS(sp_info, file = paste0(sdm_path, "fourmi_bourdon/sp_info_with_var.RDS"))



sel_vars <- readRDS(paste0(sdm_path, "fourmi_bourdon/selected_variables.RDS"))
sp_info <- readRDS(paste0(sdm_path, "fourmi_bourdon/sp_info_with_var.RDS"))


hist(sp_info$nb_variables)
table(sp_info$occ_qty, sp_info$nb_variables)

sp_info005 <- readRDS(paste0(sdm_path, "var_selection/04_sp_info_with_var_imp005.RDS"))
hist(sp_info005$nb_variables)
table(sp_info005$occ_qty, sp_info$nb_variables)

sp_info56 <- readRDS(paste0(sdm_path, "var_selection/04_sp_info_with_var_imp005_n6.RDS"))
hist(sp_info56$nb_variables)
table(sp_info56$occ_qty, sp_info$nb_variables)

