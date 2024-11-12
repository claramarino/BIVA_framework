# 1. Get the final vulnerability for all normalization methods
# 2. Compare the ouptut of each method
# Normalization methods are:
#   - max-min scaling (main text)
#   - cumulative ranking
#   - log-transformed

rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(biscale)

# select resolution
res = "110" # 55km or 110km

# load grid 
grid <- readRDS("Data/derived-data/04_Grid_110km.rds")
grid$cell_id = 1:nrow(grid)
range(grid$grid_id)

# load data
dl_expo_norm <- readRDS("Data/data-for-analyses/10_Exposure_normalized_110_km.rds")
dl_sensit <- readRDS("Data/data-for-analyses/13_Sensitivity_normalized_110_km.rds")



# What are the high/low vulnerability zones
# Are they identical based on the normalization methods?

# Combine E ans S to get Vu categories

vu_data <- function(class, norm_method){
  if(!norm_method %in% c("log","max_min", "rank")){
    stop("norm must be 'log', 'rank', or 'max_min'")}
  
  dat1 <- dl_sensit[[class]] 
  dat2 <- dl_expo_norm[[class]]
  
  if(norm_method=="max_min"){
    dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
      replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
      mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
      replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
      rename(Sensitivity = SR_ias_a_max_min)}
  
  if(norm_method=="log"){
    dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
      replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
      mutate(Exposure = expo_log/max(expo_log)) %>%
      replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
      rename(Sensitivity = SR_ias_a_log)}
  
  if(norm_method=="rank"){
    dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
      replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
      mutate(Exposure = expo_rank/max(expo_rank)) %>%
      replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
      rename(Sensitivity = SR_ias_a_rank)}
  
  data <- bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 3)
  
  data_vu <- data %>%
    mutate(vu = case_when(
      bi_class=="1-1" ~ "VL",
      bi_class %in% c("1-2", "2-1") ~ "L",
      bi_class %in% c("1-3", "3-1", "2-2") ~ "M",
      bi_class %in% c("3-2", "2-3") ~ "H",
      bi_class=="3-3" ~ "VH")) %>%
    mutate(vu = factor(vu, ordered = T, levels = c("VL","L","M","H","VH"))) %>%
    select(grid_id, cell_id, Sensitivity, Exposure, bi_class, vu) %>%
    rename(Vulnerability = vu)
  
  return(data_vu)
}

# set color code for VU categories
color_code <- c("VL" = "#FEECEB", "L" = "#FCC7C3", "M" = "#F88279", 
                "H" = "#D63B2F", "VH" = "#7A221B")


# plot the distribution of Vu categories accroding to E and S values
# map the final Vu zones for each normalization method

plot_norm_methods <- function(class){
  # get data with each normalization method
  vu_log <- vu_data(class, "log")
  table(vub_log$Vulnerability)
  vu_mm <- vu_data(class, "max_min")
  table(vub_mm$Vulnerability)
  vu_rk <- vu_data(class, "rank")
  table(vub_rk$Vulnerability)
  
  # plot distribution of VU categories
  eslog <- ggplot(vu_log)+
    geom_point(aes(x = Exposure, y = Sensitivity, color = Vulnerability))+
    scale_color_manual(values = color_code)+
    labs(title="Log-transformed")+ theme_classic()
  esmm <- ggplot(vu_mm)+
    geom_point(aes(x = Exposure, y = Sensitivity, color = Vulnerability))+
    scale_color_manual(values = color_code)+
    labs(title="Max-min scaling")+ theme_classic()
  esrk <- ggplot(vu_rk)+
    geom_point(aes(x = Exposure, y = Sensitivity, color = Vulnerability))+
    scale_color_manual(values = color_code)+
    labs(title="Cumulative ranking")+ theme_classic()
  # plot map with VU for each normalization method
  mm <- ggplot(vu_mm)+
    geom_sf(aes(fill=Vulnerability), color=NA) +
    scale_fill_manual(values = color_code)+
    theme_classic()
  rk <- ggplot(vu_rk)+
    geom_sf(aes(fill=Vulnerability), color=NA) +
    scale_fill_manual(values = color_code)+
    theme_classic()
  log <- ggplot(vu_log)+
    geom_sf(aes(fill=Vulnerability), color=NA) +
    scale_fill_manual(values = color_code)+
    theme_classic()
  
  a = ggarrange(ggarrange(esmm, esrk, eslog, ncol=3, nrow=1, legend = "none"),
                ggarrange(mm, rk, log, ncol=3, nrow=1, common.legend = T, legend = "bottom"),
                nrow=2, ncol=1)
  return(a)
}

# save the figures for supplementary materials

png("Fig/Suppl_Fig1A_Normalization_VU_bird.png", 
    height = 5 , width = 8, units = "in", res=300)
plot_norm_methods("bird")
dev.off()
png("Fig/Suppl_Fig1B_Normalization_VU_mam.png", 
    height = 5 , width = 8, units = "in", res=300)
plot_norm_methods("mam")
dev.off()
png("Fig/Suppl_Fig1C_Normalization_VU_rept.png", 
    height = 5 , width = 8, units = "in", res=300)
plot_norm_methods("rept")
dev.off()

# contingency tables for comparing normalization methods
# what are the differences in categories between max-min scaling and other methods?

contingency_plot <- function(class){
  # get data with each normalization method
  vu_log <- vu_data(class, "log")
  vu_mm <- vu_data(class, "max_min")
  vu_rk <- vu_data(class, "rank")
  
  # min max and ranking
  rk_mm <- left_join(vu_rk %>% st_drop_geometry() %>% 
                       select(grid_id, Vulnerability) %>% rename(vu_rk = Vulnerability),
                     vu_mm %>% st_drop_geometry() %>% 
                       select(grid_id, Vulnerability) %>% rename(vu_mm = Vulnerability))
  # min max and log 
  lg_mm <- left_join(vu_log %>% st_drop_geometry() %>% 
                       select(grid_id, Vulnerability) %>% rename(vu_lg = Vulnerability),
                     vu_mm %>% st_drop_geometry() %>% 
                       select(grid_id, Vulnerability) %>% rename(vu_mm = Vulnerability))
  
  lg <- ggplot(lg_mm)+
    geom_point(aes(x=vu_mm, y = vu_lg, color = vu_mm), 
               position = "jitter", alpha=.5)+
    scale_color_manual(values = color_code)+
    theme_bw()+ xlab("Max-min scaling") + ylab("Log-transformed")
  rk <- ggplot(rk_mm)+
    geom_point(aes(x=vu_mm, y = vu_rk, color = vu_mm), 
               position = "jitter", alpha=.5)+
    scale_color_manual(values = color_code)+
    theme_bw()+ xlab("Max-min scaling") + ylab("Cumulative ranking")
  
  ggarrange(rk, lg, nrow=1, ncol = 2, legend = "none")
  
}

# save the figures for supplementary materials

png("Fig/Suppl_Fig2A_Cate_Norm_bird.png", 
    height = 3, width = 6, units = "in", res=300)
contingency_plot("bird")
dev.off()
png("Fig/Suppl_Fig2B_Cate_Norm_mam.png", 
    height = 3, width = 6, units = "in", res=300)
contingency_plot("mam")
dev.off()
png("Fig/Suppl_Fig2C_Cate_Norm_rept.png", 
    height = 3, width = 6, units = "in", res=300)
contingency_plot("rept")
dev.off()

