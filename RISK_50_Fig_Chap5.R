# Figures for Chapter 5


rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(cartogram)

# select resolution
res = "110" # 55km or 110km

# load grid 
grid <- readRDS(paste0("Output/RISK_32_grid_",res, "km"))
range(grid$grid_id)

# load data
dl_expo_norm <- readRDS(paste0("Output/Exposure/RISK_23_expo_norm_", res, "_km"))
dl_sensit <- readRDS(paste0("Output/Sensitivity/RISK_33_sensit_norm_", res, "_km"))
df_compl <- readRDS(paste0("Output/Exposure/RISK_22_Comp_smpl_eff_", res, "km"))


#### Suppl Fig 1 correlation norm method ####

# for exposure
expo_norm <- list()
for(i in 1:length(dl_expo_norm)){
  expo_norm[[i]] <- dl_expo_norm[[i]] %>%
    mutate(class = names(dl_expo_norm)[i]) %>%
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

lapply(dl_expo_norm, function(x){
  cor(x %>% select(contains("expo")), method = "spearman")
})

hist(dl_expo_norm$rept$expo_max_min)
hist(dl_expo_norm$bird$expo_log)
hist(dl_expo_norm$bird$expo_rank)

hist(dl_expo_norm$mam$expo_max_min)
hist(dl_expo_norm$mam$expo_log)
hist(dl_expo_norm$mam$expo_rank)

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


#### Fig 2 Tot expo + Completion expo ####

# aggregate to larger grain size
# sinon on ne voit pas la distorsion

df_5deg <- grid %>% 
  st_make_grid(cellsize = 550000, square = F) %>% 
  st_sf()
ggplot(df_5deg)+geom_sf()

df_compl_5deg <- st_intersection(df_compl, 
                                df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
  st_drop_geometry() %>%
  group_by(id_5deg) %>%
  summarise(comp_sum = mean(na.omit(comp_se_sum)),
            comp_prod = mean(na.omit(comp_se_product)))

df_compl_5deg_sf <- left_join(df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg)),
                             df_compl_5deg)

ggplot(df_compl_5deg_sf)+geom_sf(aes(fill=comp_prod))


cartog_prod <- cartogram_ncont(df_compl_5deg_sf,
                               weight = "comp_prod")

# first look cartogram
ggplot(cartog_prod)+
  geom_sf(aes(fill=comp_prod), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)

# add exposure value for color
expo <- left_join(
  grid, 
  dl_expo_norm$all_groups %>% dplyr::select(grid_id, SR_tot_ias_log)) %>%
  rename(SR_tot_ias = SR_tot_ias_log) %>%
  replace_na(list(SR_tot_ias=0))

# check color
ggplot(expo)+
  geom_sf(aes(fill=SR_tot_ias), color = NA)

# aggregate expo at higher resolution

df_expo_5deg <- st_intersection(expo,
                                df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
  st_drop_geometry() %>%
  group_by(id_5deg) %>%
  summarise(SR_ias = mean(na.omit(SR_tot_ias)))


cartog_compl_expo <- left_join(cartog_prod, df_expo_5deg)

world <- rnaturalearth::ne_countries(returnclass = "sf")
world_moll <- st_make_valid(world)%>%
  st_transform(crs = st_crs(grid))

ggplot()+
  geom_sf(data = grid, fill = "grey40", color = NA)+
  geom_sf(data = cartog_compl_expo, aes(fill=SR_ias), color = NA) +
  scale_fill_viridis_c(option="magma", direction = -1)+
  theme_classic()


##### Fig 3 Sensitivity ~ Total SR ####

sr_all <- readRDS(paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_ias_a_t_nt_", res))

grid_terr$cell_id = 1:nrow(grid_terr)

# mammals 

grid_mam <- left_join(grid_terr, sr_all$mam) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

lm_mam <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_mam)
grid_mam$res = residuals(lm_mam)

m<- ggplot(grid_mam)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Mammalia")

# birds 

grid_bird <- left_join(grid_terr, sr_all$bird) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))
lm_bird <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_bird)
grid_bird$res = residuals(lm_bird)

b<-ggplot(grid_bird)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Aves")

# reptiles

grid_rept <- left_join(grid_terr, sr_all$rept) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

lm_rept <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_rept)
grid_rept$res = residuals(lm_rept)

r <-  ggplot(grid_rept)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Reptilia")

pdf(paste0("Fig/RISK_33_resid_iasA_SRtot_BMR_", res, ".pdf"), 8, 5)
ggarrange(b, m, r, ncol = 2, nrow=2, common.legend = T)
dev.off()
