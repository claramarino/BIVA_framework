# Figures for Chapter 5


rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)


# select resolution
res = "55" # 55km or 110km

# load grid 
grid <- readRDS(paste0("Output/RISK_32_grid_",res, "km"))
range(grid$grid_id)

# load data
dl_all_norm <- readRDS(paste0("Output/Exposure/RISK_23_expo_norm_", res, "_km"))

#### Suppl Fig 1 correlation norm method ####

# for exposure
expo_norm <- list()
for(i in 1:length(dl_all_norm)){
  expo_norm[[i]] <- dl_all_norm[[i]] %>%
    mutate(class = names(dl_all_norm)[i]) %>%
    select(grid_id, class, contains("expo"))
}

df_expo_norm <- bind_rows(expo_norm)

log_max_min <- ggplot(df_expo_norm %>% filter(class!="all_groups"),
       aes(x=expo_max_min, y = expo_log)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'lm')+
  theme_classic()


log_rank <- ggplot(df_expo_norm %>% filter(class!="all_groups"),
                      aes(x=expo_rank, y = expo_log))+
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = "lm")+
  theme_classic()


rank_max_min <- ggplot(df_expo_norm %>% filter(class!="all_groups"),
                      aes(x=expo_max_min, y = expo_rank)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'lm')+
  theme_classic()


ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)

lapply(dl_all_norm, function(x){
  cor(x %>% select(contains("expo")), method = "spearman")
})

hist(dl_all_norm$rept$expo_max_min)
hist(dl_all_norm$bird$expo_log)
hist(dl_all_norm$bird$expo_rank)

hist(dl_all_norm$mam$expo_max_min)
hist(dl_all_norm$mam$expo_log)
hist(dl_all_norm$mam$expo_rank)

# on prend expo max-min => la + simple
# faudra juste vérifier que les bivariate sont identiques si on prend log

# for sensitivity



#### Suppl Fig 2 correlation metrics expo ####

# select normalization method
norm = "rank"# max_min, log or rank
metrics <- list()
for(i in 1:length(dl_all_norm)){
  metrics[[i]] <- dl_all_norm[[i]] %>%
    mutate(class = names(dl_all_norm)[i]) %>%
    select(grid_id, class, contains(norm))
}

df_metrics <- bind_rows(metrics)

colnames(df_metrics) <- c(
  "grid_id",  "class", "SR_tot_ias", 
  "med_range","med_ib", "expo" )

sr_range <- ggplot(df_metrics %>% filter(class!="all_groups"),
                      aes(x=SR_tot_ias, y = med_range))  +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'lm')+
  theme_classic()

imp_range <- ggplot(df_metrics %>% filter(class!="all_groups"),
                   aes(x=med_ib, y = med_range)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'lm')+
  theme_classic()


imp_sr <- ggplot(df_metrics %>% filter(class!="all_groups"),
                       aes(x=med_ib, y = SR_tot_ias))  +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'lm')+
  theme_classic()


ggarrange(sr_range, imp_range, imp_sr,
          ncol = 2, nrow = 2, common.legend = T)

lapply(dl_metrics, function(x){
  cor(x %>% select(-grid_id), 
      method = "pearson")
})





#### Fig 1 Expo + sensib = risk #####  



#### Fig 2 Tot expo + Bias expo ####





##### ####