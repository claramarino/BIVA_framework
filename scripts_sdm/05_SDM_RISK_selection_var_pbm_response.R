# pbm selection var

rm(list = ls())

#library(raster)
#devtools::install_github("biomodhub/biomod2", dependencies = TRUE)
library(biomod2)
library(tidyverse)
library(Rarity)
library(terra)

# file path
sdm_path <- "Z:/R/sdmrisk/"


sp_info <- readRDS(paste0(sdm_path, "var_selection/04_sp_info_with_var.RDS"))


sp1 = "sp1314773" # fourmi
sp2 = "sp1340503" # bourdon

sp=sp1

sp_info %>% filter(sp_key == substring(sp, 3))

model_runs <- readRDS(paste0(sdm_path, "var_selection/all_vars/", sp, "/model_runs.RDS"))
input_data <- get(load(model_runs@formated.input.data@link))

# Calcul des courbes de réponse
resp <- bm_PlotResponseCurves(bm.out = model_runs,
                              fixed.var = "mean",
                              data_species = input_data@data.species
)$tab

colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")

p <- ggplot(resp, aes(x = Var.value, y = Response))+ 
  geom_line(alpha = 0.2, aes(group = Model)) + 
  stat_smooth() +
  facet_wrap(~Variable, scales = "free_x") + 
  theme_bw() + 
  ylim(0, 1.1) + 
  xlab("Variable value")

respf <- resp %>% filter(Variable %in% c("pop","bio16","bio5","bio6","bio4")) %>%
  mutate(Variable = factor(Variable, levels = c("pop","bio16","bio5","bio6","bio4"),
                           ordered = T))

ggplot(respf %>% filter(Variable %in% c("bio5", "bio6", "bio4")), aes(x = Var.value, y = Response))+ 
  geom_line(alpha = 0.8, aes(group = Model, color = Model))+
  facet_wrap(~Variable, scales = "free_x") + 
  theme_bw() + 
  ylim(0, 1.1) + 
  xlab("Variable value")


# bourdon 

sp = sp2
model_runs <- readRDS(paste0(sdm_path, "var_selection/all_vars/", sp, "/model_runs.RDS"))
input_data <- get(load(model_runs@formated.input.data@link))

# Calcul des courbes de réponse
resp <- bm_PlotResponseCurves(bm.out = model_runs,
                              fixed.var = "mean",
                              data_species = input_data@data.species
)$tab

colnames(resp) <- c("Index", "Variable", "Var.value", "Model", "Response")





respf <- resp %>% filter(Variable %in% c("bio6","bio8","bio5","bio4")) %>%
  mutate(Variable = factor(Variable, levels = c("bio6","bio8","bio5","bio4"),
                           ordered = T))
ggplot(respf, aes(x = Var.value, y = Response))+ 
  geom_line(aes(group = Model, color = Model), linewidth=1)+
  facet_wrap(~Variable, scales = "free_x") + 
  theme_bw() + 
  ylim(0, 1.1) + 
  xlab("Variable value")


pres_bg <- bind_cols(data.frame(pres = input_data@data.species),
                     input_data@data.env.var)

df <- pres_bg %>% select(pres, bio6, bio8, bio5, bio4) %>%
  mutate(pres = as.character(pres)) %>%
  pivot_longer(cols = bio6:bio4,
               names_to = "var", 
               values_to = "value")

ggplot(df)+
  geom_boxplot(aes(x=pres, y = value))+
  facet_wrap(~var, scales = "free")

# overlap 
coord <- bind_cols(data.frame(pres = as.character(input_data@data.species)),
                   input_data@coord)

ggplot(coord)+
  geom_point(aes(x=x, y=y, color = pres), alpha = .3)

# A réfléchir : ajouter limites des données de présence ?

png(paste0("graphiques/response_plot_", sp, ".png"), width = 550 * 4.2, height = 550 * 4.2, res = 300)

print(p)

dev.off()