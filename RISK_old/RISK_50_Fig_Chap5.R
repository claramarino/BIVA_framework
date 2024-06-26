# Figures for Chapter 5


rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(cartogram)
library(cowplot)
library(grid)
library(gridExtra) 
library(biscale)

# select resolution
res = "110" # 55km or 110km

# load grid 
grid <- readRDS(paste0("Output/RISK_32_grid_",res, "km"))
grid$cell_id = 1:nrow(grid)
range(grid$grid_id)

# load data
dl_expo_norm <- readRDS(paste0("Output/Exposure/RISK_23_expo_norm_", res, "_km"))
dl_sensit <- readRDS(paste0("Output/Sensitivity/RISK_33_sensit_norm_", res, "_km"))
df_compl <- readRDS(paste0("Output/Exposure/RISK_22_Comp_smpl_eff_", res, "km"))
df_compl_sens <- readRDS(paste0("Output/Sensitivity/Completeness/RISK_34_completeness_ed_se_", res))

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
  geom_smooth(aes(color = class), method = 'loess')+
  theme_classic()


log_rank <- ggplot(df_expo_norm %>% filter(class!="all_groups"),
                      aes(x=expo_rank, y = expo_log))+
  geom_point(color = "grey", alpha = .2)+
  geom_smooth(aes(color = class), method = "loess")+
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

pdf("Fig/Suppl_Fig3A_Norm_exposure.pdf", 6, 5)
ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)
dev.off()


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


pdf("Fig/Suppl_Fig3B_Norm_sensitivity.pdf", 6, 5)
ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)
dev.off()

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

#### Suppl Fig 3 spatial diff metrics expo ####
grid_all <- left_join(grid, dl_expo_norm$all)


all_sr <- ggplot(grid_all)+
  geom_sf(aes(fill = SR_tot_ias_max_min), color = NA)+
  scale_fill_viridis_c(option="viridis", direction = -1, na.value = "grey80" )+
  theme_classic()
all_impact <- ggplot(grid_all)+
  geom_sf(aes(fill = med_ib_max_min), color = NA)+
  scale_fill_viridis_c(option="viridis", direction = -1, na.value = "grey80" )+
  theme_classic()
all_range <- ggplot(grid_all)+
  geom_sf(aes(fill = med_range_max_min), color = NA)+
  scale_fill_viridis_c(option="viridis", direction = -1, na.value = "grey80" )+
  theme_classic()

pdf("Fig/Chap5_FigS3_spatial_expo_metrics.pdf", 8, 12)
ggarrange(all_sr, all_range, all_impact, ncol = 1, nrow=3,
          common.legend = T, legend = "top")
dev.off()
#### Fig 1 Expo + sensib = risk #####  

# cb de cells avec au moins une des 304 IAS
lapply(dl_expo_norm, function(x){
  nrow(x)/17257
})


# cb de cells avec sp sensibles
lapply(df_compl_sens, function(x){
  nrow(x %>% filter(SR_ias_a>0))/17257
})


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
                   style = "fisher", dim = 3)
  
  # bi_class_breaks(dat_sf, x = Exposure, y = Sensitivity, style = "fisher",
  #                 dim = 4, dig_lab = c(4, 5), split = FALSE)
  
  map <- ggplot() +
    geom_sf(data = data, aes(fill = bi_class), 
            color = NA, show.legend = FALSE) + #
    bi_scale_fill(pal = "BlueOr", dim = 3)+
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
biv
pdf("Fig/Chap5_Fig1_expo_sensib_biv.pdf", 8,6)
ggarrange(e, s, biv, ncol = 3, common.legend = F)
dev.off()


legend_biv <- bi_legend(pal = "BlueOr",
                    dim = 3,
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

# correlation ed se 
df_com <- df_compl %>% st_drop_geometry() %>%
  filter(!is.na(norm_dens))
cor.test(df_com$norm_dens, df_com$DB_and_GRIIS)

ggplot(df_com, aes(x=norm_dens, y =DB_and_GRIIS))+
  geom_point()+
  geom_smooth(method="lm")

mod <- lm(norm_dens~DB_and_GRIIS, data= df_com)
summary(mod)
plot(mod)


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

mod = lm(expo_norm~comp_se_product+class, data = rel_expo_comp)
summary(mod)


mod = lm(SR_tot_ias_log~comp_se_product+class, data = rel_expo_comp)
summary(mod)


plot(mod)
cor.test(rel_expo_comp$SR_tot_ias_log, rel_expo_comp$comp_se_product)
cor.test(rel_expo_comp$expo_norm, rel_expo_comp$comp_se_product)


#### Fig 4 Completion sensitivity ####


# range de valeur completeness
lapply(df_compl_sens, summary)

grid_compl_sens <- left_join(grid, df_compl_sens$rept)

# relationship between both components of completeness
ggplot(df_compl_sens$all, 
       aes(x = se, y = ed))+
  geom_point(alpha = .3, color = "grey60", size=2)+
  xlab("Sampling effort IUCN")+
  ylab("Knowledge on IAS threat")+
  theme_classic()

# plot iucn sampling effort 
ggplot(grid_compl_sens)+
  geom_sf(aes(fill=se), color=NA) +
  scale_fill_viridis_c(na.value = "grey80")+
  theme_classic()+
  ggtitle("IUCN sampling effort")

# global knowledge ias threat
ggplot(grid_compl_sens)+
  geom_sf(aes(fill=ed), color=NA) +
  scale_fill_viridis_c(na.value = "grey80")+
  theme_classic()+
  ggtitle("Knowledge on IAS threat")

#final completeness = same as ed
ggplot(grid_compl_sens)+
  geom_sf(aes(fill=comp_prod), color=NA) +
  scale_fill_viridis_c(na.value = "grey80")+
  theme_classic()+
  ggtitle("Completeness")



########### try cartogram sensit with completeness ?

# cartogram function for each group
class="rept"

# initialize grid 5deg
df_5deg <- grid %>% 
  st_make_grid(cellsize = 550000, square = F) %>% 
  st_sf()

cartog_sensit <- function (class){
  # get grid completeness for the class
  grid_compl_sens <- left_join(grid, df_compl_sens[[class]])
  
  #intersection with 5deg grid
  df_compl_5deg <- st_intersection(grid_compl_sens, 
                                   df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
    st_drop_geometry() %>%
    group_by(id_5deg) %>%
    summarise(comp_prod = mean(na.omit(comp_prod)))
  df_compl_5deg_sf <- left_join(df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg)),
                                df_compl_5deg)
  
  cartog_prod <- cartogram_ncont(df_compl_5deg_sf,
                                 weight = "comp_prod")
  
  # add sensitivity value for color
  sensit <- left_join(
    grid,
    dl_sensit[[class]] %>% dplyr::select(cell_id, SR_ias_a_max_min)) %>%
    rename(sensitivity = SR_ias_a_max_min) %>%
    replace_na(list(sensitivity=0))
  # aggregate sensit at higher resolution
  df_sensit_5deg <- st_intersection(sensit,
                                  df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
    st_drop_geometry() %>%
    group_by(id_5deg) %>%
    summarise(Sensitivity = mean(na.omit(sensitivity)))
  
  
  cartog_compl_expo <- left_join(cartog_prod, df_sensit_5deg)
  
  
  ggplot()+
    geom_sf(data = grid, fill = "grey80", color = NA)+
    geom_sf(data = cartog_compl_expo, aes(fill=Sensitivity), color = NA) +
    scale_fill_viridis_c(option="viridis")+
    theme_classic()
}


cb <- cartog_sensit("bird")
cm <- cartog_sensit("mam")
cr <- cartog_sensit("rept")


cb
cm
cr
pdf("Fig/Fig4_Sensit_compl_map_BMR.pdf", 9,8)
ggarrange(cb, cm, cr, legend = "top", ncol = 2, nrow=2)
dev.off()



# sensitivity ~ completeness 


sens_compl <- function(class){
  compl_sens <- left_join(
    dl_sensit[[class]] %>% dplyr::select(cell_id, SR_ias_a_max_min),
    df_compl_sens[[class]]) %>%
  rename(sensitivity = SR_ias_a_max_min) %>%
    filter(sensitivity>0)

  # ggplot(compl_sens, 
  #        aes(x = comp_prod, y = sensitivity))+
  #   geom_point(alpha = .3, color = "grey60", size=2)+
  #   #geom_smooth(aes(color=Class), method = "lm")+
  #   xlab("Completeness")+
  #   ylab("Sensitivity")+
  #   geom_rug(size=0.1) +
  #   theme_classic()+
  #   theme(legend.position = "top")
  ggscatterhist(
    compl_sens, x = "comp_prod", y = "sensitivity",
    color = "maroon", alpha = .3,
    margin.plot = "density",
    margin.params = list(fill = "lightpink", color = "maroon"),
    xlab = "Completeness",
    ylab = "Sensitivity",
    main.plot.size =1,
    ggtheme = theme_classic()
  )
    
}


pdf("Fig/Chap5_Fig4_Sensit_compl_bird.pdf", 3,3)
sens_compl("bird")
dev.off()

pdf("Fig/Chap5_Fig4_Sensit_compl_mam.pdf", 3,3)
sens_compl("mam")
dev.off()

pdf("Fig/Chap5_Fig4_Sensit_compl_rept.pdf", 3,3)
sens_compl("rept")
dev.off()


ggarrange(scb, scm, scr, ncol = 3, nrow=1)

##### Fig 3 Sensitivity ~ Total SR ####

sr_all <- readRDS(paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_ias_a_t_nt_", res))

grid$cell_id = 1:nrow(grid)

# mammals 

grid_mam <- left_join(grid, sr_all$mam) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

lm_mam <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_mam)
summary(lm_mam)

cor.test(log(grid_mam$SR_ias_a+1), log(grid_mam$SR_tot+1))


grid_mam$res = residuals(lm_mam)

m<- ggplot(grid_mam)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Mammalia")

# birds 

grid_bird <- left_join(grid, sr_all$bird) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))
lm_bird <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_bird)
grid_bird$res = residuals(lm_bird)

cor.test(log(grid_bird$SR_ias_a+1), log(grid_bird$SR_tot+1))

b<-ggplot(grid_bird)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Aves")

# reptiles

grid_rept <- left_join(grid, sr_all$rept) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

lm_rept <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_rept)
grid_rept$res = residuals(lm_rept)

cor.test(log(grid_rept$SR_ias_a+1), log(grid_rept$SR_tot+1))

r <-  ggplot(grid_rept)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Reptilia")

pdf(paste0("Fig/RISK_33_resid_iasA_SRtot_BMR_", res, ".pdf"), 8, 5)
ggarrange(b, m, r, ncol = 2, nrow=2, common.legend = T)
dev.off()


###### Fig 4 - High risk cells only #######

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


rb <- plot_high_risk("bird")+
  labs(title = "Aves")
rm <- plot_high_risk("mam")+
  labs(title = "Mammalia")
rr <- plot_high_risk("rept")+
  labs(title = "Reptilia")

pdf("Fig/Chap5_Suppl_fg_high_risk.pdf", 6, 5)
ggarrange(rb, rm, rr, ncol= 2, nrow=2, common.legend = T)
dev.off()

# plot high risk and low risk zones

bi_pal(pal = "BlueOr",
       dim = 4,preview = F)
class ='bird'
plot_highlow_risk <- function(class){
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
    mutate(high_risk = if_else(bi_class=="3-3", "1", "other")) %>%
    mutate(high_risk = if_else(bi_class=="1-1", "0", high_risk))
  
  table(data_high$high_risk)
  ggplot(data_high)+
    geom_sf(aes(fill=high_risk), color=NA) +
    scale_fill_manual(values = c("1" = "green4", "0" = "wheat1", "other" = "palegreen2"))+
    theme_classic()
}


rb <- plot_highlow_risk("bird")+
  labs(title = "Aves")
rm <- plot_highlow_risk("mam")+
  labs(title = "Mammalia")
rr <- plot_highlow_risk("rept")+
  labs(title = "Reptilia")

pdf("Fig/Chap5_Fig2_high_and_low_risk.pdf", 6, 5)
ggarrange(rb, rm, rr, ncol= 2, nrow=2, common.legend = T)
dev.off()



ggplot(data_high)+
  geom_sf(aes(fill=high_risk), color=NA) +
  scale_fill_manual(values = c("1" = "maroon", "0" = "wheat1", "other" = "lightpink"))+
  theme_classic()


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

pdf("Fig/Fig4_vulnerability_tetrapods.pdf", 6, 5)
ggarrange(vb, vm, vr, ncol= 2, nrow=2, common.legend = T)
dev.off()




# VU as the combination of E and S

e <- ggarrange(eb, em, er,
               ncol = 1, nrow = 3, legend = "none")
s <- ggarrange(sb, sm, sr,
               ncol = 1, nrow = 3, legend = "none")
vu <- ggarrange(vb, vm, vr, 
                ncol= 1, nrow=3, legend = "none")

pdf("Fig/Fig2_Expo_sensib_VU.pdf", 8,6)
ggarrange(e, s, vu, ncol = 3, common.legend = F)
dev.off()

legend_vu <- cowplot::get_legend(vr + theme(legend.direction = "horizontal"))


pdf("Fig/Fig2_Legend_vu.pdf")
grid.newpage()
grid.draw(legend_vu)
dev.off()


