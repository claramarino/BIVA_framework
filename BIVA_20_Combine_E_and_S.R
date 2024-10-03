# Combine Exposure and Sensitivity to get final vulnerability

rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(biscale)

# select resolution
res = "110" # 55km or 110km

# load grid 
grid <- readRDS(paste0("Output/RISK_32_grid_",res, "km"))
grid$cell_id = 1:nrow(grid)
range(grid$grid_id)

# load data
dl_expo_norm <- readRDS(paste0("Output/Exposure/BIVA_10_expo_norm_", res, "_km"))
dl_sensit <- readRDS(paste0("Output/Sensitivity/BIVA_13_sensit_norm_", res, "_km"))
df_compl <- readRDS(paste0("Output/Exposure/RISK_22_Comp_smpl_eff_", res, "km"))
df_compl_sens <- readRDS(paste0("Output/Sensitivity/Completeness/RISK_34_completeness_ed_se_", res))

# openxlsx::write.xlsx(df_compl %>% sf::st_drop_geometry(), 
#                      file = paste0("Output/Exposure/BIVA_20_expo_completeness_", res, "_km.xlsx"))
# 
# openxlsx::write.xlsx(df_compl_sens, 
#                      file = paste0("Output/Sensitivity/BIVA_20_sensi_completeness_", res, "_km.xlsx"))


ntot = nrow(grid)

# percentage of cells with at least one of the 304 IAS
lapply(dl_expo_norm, function(x){
  nrow(x)/ntot
})

# percentage of cells with at least one sensitive native species to the 304 IAs
lapply(df_compl_sens, function(x){
  nrow(x %>% filter(SR_ias_a>0))/17257
})


plot_expo <- function(class){
  # class is bird, mam, rept
  if(!class %in% c("bird","mam","rept")){
    stop("class must be 'bird', 'mam', or 'rept")}
  dat <- dl_expo_norm[[class]] 
  dat_sf <- left_join(grid, dat) %>%
    select(grid_id, cell_id, contains("expo")) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min))
  p <- ggplot(dat_sf)+
    geom_sf(aes(fill = Exposure), color = NA) +
    scale_fill_viridis_c(option="inferno", direction = -1)+
    theme_classic()
  return(p)
}

eb <- plot_expo("bird")
em <- plot_expo("mam")
er <- plot_expo("rept")
legend_e <- cowplot::get_legend(er+ theme(legend.direction = "horizontal"))


plot_sensib <- function(class){
  dat <- dl_sensit[[class]] 
  dat_sf <- left_join(grid, dat) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  p <- ggplot(dat_sf)+
    geom_sf(aes(fill = Sensitivity), color = NA) +
    scale_fill_viridis_c(option="mako", direction = -1)+
    theme_classic()
  return(p)
}

sb <- plot_sensib("bird")
sm <- plot_sensib("mam")
sr <- plot_sensib("rept")
legend_s <- cowplot::get_legend(sr + theme(legend.direction = "horizontal"))


# Define low to high exposure and sensitivity
# merge the two information
# get vulnerability categories: VL, L, M, H, VH
plot_vulnerability <- function(class){
  dat1 <- dl_sensit[[class]] 
  dat2 <- dl_expo_norm[[class]]
  
  dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  
  data <- bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 3)
  
  data_vu <- data %>%
    mutate(vu = case_when(
      bi_class=="1-1" ~ "VL",
      bi_class %in% c("1-2", "2-1") ~ "L",
      bi_class %in% c("1-3", "3-1", "2-2") ~ "M",
      bi_class %in% c("3-2", "2-3") ~ "H",
      bi_class=="3-3" ~ "VH")) %>%
    mutate(vu = factor(vu, ordered = T, levels = c("VL","L","M","H","VH")))
  
  table(data_vu$vu)
  ggplot(data_vu)+
    geom_sf(aes(fill=vu), color=NA) +
    scale_fill_manual(values = c("VL" = "#FEECEB", 
                                 "L" = "#FCC7C3", 
                                 "M" = "#F88279", 
                                 "H" = "#D63B2F", 
                                 "VH" = "#7A221B"))+
    theme_classic()
}

vb <- plot_vulnerability("bird")+
  labs(title = "Aves")
vm <- plot_vulnerability("mam")+
  labs(title = "Mammalia")
vr <- plot_vulnerability("rept")+
  labs(title = "Reptilia")


# VU as the combination of E and S

e <- ggarrange(eb, em, er,
               ncol = 1, nrow = 3, legend = "none")
s <- ggarrange(sb, sm, sr,
               ncol = 1, nrow = 3, legend = "none")
vu <- ggarrange(vb, vm, vr, 
                ncol= 1, nrow=3, legend = "none")

ggarrange(e, s, vu, ncol = 3, common.legend = F)


# Save final figure
pdf("Fig/Fig2_Expo_sensib_VU.pdf", 8,6)
ggarrange(e, s, vu, ncol = 3, common.legend = F)
dev.off()

#save legend for VU color code
legend_vu <- cowplot::get_legend(vr + theme(legend.direction = "horizontal"))
pdf("Fig/Fig2_Legend_vu.pdf")
grid.newpage()
grid.draw(legend_vu)
dev.off()
