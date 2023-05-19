# Figures for Chapter 5


rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(cartogram)
library(cowplot)
library(grid)
library(gridExtra) 

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
    select(grid_id, class, contains("expo"), SR_tot_ias_log)
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
  round(cor(x %>% select(contains("expo")), method = "spearman"), 3)
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


sensi_norm <- list()
for(i in 1:length(dl_sensit)){
  sensi_norm[[i]] <- dl_sensit[[i]] %>%
    mutate(class = names(dl_sensit)[i])
}

df_sensi_norm <- bind_rows(sensi_norm)


log_max_min <- ggplot(df_sensi_norm,
                      aes(x=SR_ias_a_max_min, y = SR_ias_a_log)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class))+
  theme_classic()


log_rank <- ggplot(df_sensi_norm,
                   aes(x=SR_ias_a_rank, y = SR_ias_a_log))+
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class))+
  theme_classic()
#log_rank

rank_max_min <- ggplot(df_sensi_norm,
                       aes(x=SR_ias_a_max_min, y = SR_ias_a_rank)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class))+
  theme_classic()


ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)


lapply(dl_sensit, function(x){
  round(cor(x %>% select(contains("SR")), method = "pearson"), 3)
})

#### Suppl Fig 2 correlation metrics expo ####

# select normalization method
norm = "max_min"# max_min, log or rank
metrics <- list()
for(i in 1:length(dl_expo_norm)){
  metrics[[i]] <- dl_expo_norm[[i]] %>%
    mutate(class = names(dl_expo_norm)[i]) %>%
    select(grid_id, class, contains(norm))
}

df_metrics <- bind_rows(metrics) %>%
  mutate(class = case_when(
    class=="bird" ~ "Aves",
    class=="mam" ~ "Mammalia",
    class=="rept" ~ "Reptilia"
  ))

colnames(df_metrics) <- c(
  "grid_id",  "class", "SR_tot_ias", 
  "med_range","med_ib", "expo" )

sr_range <- ggplot(df_metrics %>% filter(class!="all_groups"),
                      aes(x=SR_tot_ias, y = med_range))  +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'loess')+
  ylab("Median alien range") +xlab("Alien species richness")+
  theme_classic()

imp_range <- ggplot(df_metrics %>% filter(class!="all_groups"),
                   aes(x=med_ib, y = med_range)) +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'loess')+
  xlab("Median impact breadth") +ylab("Median alien range")+
  theme_classic()


imp_sr <- ggplot(df_metrics %>% filter(class!="all_groups"),
                       aes(x=med_ib, y = SR_tot_ias))  +
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = 'loess')+
  xlab("Median impact breadth") +ylab("Alien species richness")+
  theme_classic()

pdf("Fig/Chap5_Fig0_Rel_exposure_components.pdf", 6, 5)
ggarrange(sr_range, imp_range, imp_sr,
          ncol = 2, nrow = 2, common.legend = T)
dev.off()


lapply(dl_metrics, function(x){
  cor(x %>% select(-grid_id), 
      method = "pearson")
})



#### Fig 1 Expo + sensib = risk #####  

grid$cell_id = 1:nrow(grid)
class="bird"
plot_expo <- function(class){
  # class is bird, mam, rept
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



plot_biv  <- function(class){
  dat1 <- dl_sensit[[class]] 
  dat2 <- dl_expo_norm[[class]]
  
  dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  
  data <- bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 4)
  
  # bi_class_breaks(dat_sf, x = Exposure, y = Sensitivity, style = "fisher",
  #                 dim = 4, dig_lab = c(4, 5), split = FALSE)
  
  map <- ggplot() +
    geom_sf(data = data, aes(fill = bi_class), 
            color = NA, show.legend = FALSE) + #
    bi_scale_fill(pal = "BlueOr", dim = 4)+
    theme_classic()

  return(map)
} 

bivb <- plot_biv("bird")
bivm <- plot_biv("mam")
bivr <- plot_biv("rept")

ggarrange(eb, sb, em, sm, er, sr,
          ncol = 2, nrow = 3, common.legend = T)

e <- ggarrange(eb, em, er,
          ncol = 1, nrow = 3, legend = "none")
s <- ggarrange(sb, sm, sr,
          ncol = 1, nrow = 3, legend = "none")


biv <- ggarrange(bivb, bivm, bivr, ncol = 1, nrow= 3, legend="none")

pdf("Fig/Chap5_Fig1_expo_sensib_biv.pdf", 8,6)
ggarrange(e, s, biv, ncol = 3, common.legend = F)
dev.off()


legend_biv <- bi_legend(pal = "BlueOr",
                    dim = 4,
                    xlab = "Higher exposure ",
                    ylab = "Higher sensitivity ",
                    size = 8)
legend_biv

pdf("Fig/Chap5_Fig1_Legend_S.pdf")
grid.newpage()
grid.draw(legend_s)
dev.off()

pdf("Fig/Chap5_Fig1_Legend_E.pdf")
grid.newpage()
grid.draw(legend_e)
dev.off()

pdf("Fig/Chap5_Fig1_Legend_biv.pdf")
legend_biv
dev.off()


plot_sensi_expo <- function(class){
  dat1 <- dl_sensit[[class]] 
  dat2 <- dl_expo_norm[[class]]
  
  dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  
  data <- bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 3)
  
  # bi_class_breaks(dat_sf, x = Exposure, y = Sensitivity, style = "fisher",
  #                 dim = 4, dig_lab = c(4, 5), split = FALSE)
  
  bip <- ggplot(data = data, aes(x = Exposure, y = Sensitivity)) +
    geom_point(aes(color = bi_class), show.legend = FALSE) + #
    bi_scale_color(pal = "BlueOr", dim = 3)+
    theme_classic()
  
  return(bip)
}


seb <- plot_sensi_expo("bird")+
  labs(title = "Aves")
sem <- plot_sensi_expo("mam")+
  labs(title = "Mammalia")
ser <- plot_sensi_expo("rept")+
  labs(title = "Reptilia")

pdf("Fig/Chap5_Suppl_fg_biva_scale.pdf", 8, 3)
ggarrange(seb, sem, ser, ncol=3, nrow=1)
dev.off()


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

# ggplot(df_compl_5deg_sf)+geom_sf(aes(fill=comp_prod))

cartog_prod <- cartogram_ncont(df_compl_5deg_sf,
                               weight = "comp_prod")

# first look cartogram
# ggplot(cartog_prod)+
#   geom_sf(aes(fill=comp_prod), color = NA) +
#   scale_fill_viridis_c(option="magma", direction = -1)

# add exposure value for color
expo <- left_join(
  grid, 
  dl_expo_norm$all_groups %>% dplyr::select(grid_id, expo_max_min)) %>% #SR_tot_ias_log
  #rename(SR_tot_ias = SR_tot_ias_log) %>%
  rename(exposure = expo_max_min) %>%
  replace_na(list(exposure=0))

# check color
# ggplot(expo)+
#   geom_sf(aes(fill=SR_tot_ias), color = NA)

# aggregate expo at higher resolution
df_expo_5deg <- st_intersection(expo,
                                df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
  st_drop_geometry() %>%
  group_by(id_5deg) %>%
  summarise(Exposure = mean(na.omit(exposure)))


cartog_compl_expo <- left_join(cartog_prod, df_expo_5deg)

pdf("Fig/Chap5_Fig2_Completeness_cartogram.pdf")
ggplot()+
  geom_sf(data = grid, fill = "grey60", color = NA)+
  geom_sf(data = cartog_compl_expo, aes(fill=Exposure), color = NA) +
  scale_fill_viridis_c(option="viridis")+
  theme_classic()
dev.off()


# add relationship between exposure and completeness
rel_expo_comp <- left_join(df_expo_norm, df_compl %>% st_drop_geometry()) %>%
  filter(class!="all_groups") %>%
  mutate(Class=case_when(
    class=="bird" ~ "Aves",
    class=="mam" ~ "Mammalia",
    class=="rept" ~ "Reptilia"
  )) %>%
  mutate(expo_norm = expo_max_min/max(expo_max_min))

pdf("Fig/Chap5_Fig2b_Expo_Complet.pdf", 3, 4)
ggplot(rel_expo_comp, 
       aes(x = comp_se_product, y = expo_norm))+
  geom_point(alpha = .3, color = "grey60", size=2)+
  geom_smooth(aes(color=Class), method = "lm")+
  xlab("Completeness")+
  ylab("Exposure")+
  theme_classic()+
  theme(legend.position = "top")
dev.off()

mod = lm(SR_tot_ias_log~comp_se_product, data = rel_expo_comp)
summary(mod)

plot(mod)


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


###### Fig 4 - #######

# plot the high risk cells only

plot_high_risk <- function(class){
  dat1 <- dl_sensit[[class]] 
  dat2 <- dl_expo_norm[[class]]
  
  dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  
  data <- bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 3)

  data_high <- data %>%
    mutate(high_risk = if_else(bi_class=="3-3", "1", "0"))
  
  ggplot(data_high)+
    geom_sf(aes(fill=high_risk), color=NA) +
    scale_fill_manual(values = c("1" = "#D73029", "0" = "grey80" ))+
    theme_classic()
  }


plot_high_risk("bird")
plot_high_risk("mam")
plot_high_risk("rept")
