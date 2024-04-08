# Map completeness and explore relationship with metrics
# for exposure and sensitivity

#### Load data ####
rm(list=ls())

library(tidyverse)
library(sf)
library(ggpubr)
library(cartogram)
library(cowplot)

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



#### Exposure ####

# relationship between completeness components 
# correlation ed and se 
df_com <- df_compl %>% st_drop_geometry() %>%
  filter(!is.na(norm_dens))
cor.test(df_com$norm_dens, df_com$DB_and_GRIIS)

ggplot(df_com, aes(x=norm_dens, y =DB_and_GRIIS))+
  geom_point()+
  geom_smooth(method="lm")

mod <- lm(norm_dens~DB_and_GRIIS, data= df_com)
summary(mod)
plot(mod)

# aggregate to larger grain size (for visualizing the distorsion)
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


# add exposure value for color
expo <- left_join(
  grid, 
  dl_expo_norm$all_groups %>% dplyr::select(grid_id, expo_max_min)) %>% #SR_tot_ias_log
  #rename(SR_tot_ias = SR_tot_ias_log) %>%
  rename(exposure = expo_max_min) %>%
  replace_na(list(exposure=0))

# aggregate expo at higher resolution
df_expo_5deg <- st_intersection(expo,
                                df_5deg %>% mutate(id_5deg = 1:nrow(df_5deg))) %>%
  st_drop_geometry() %>%
  group_by(id_5deg) %>%
  summarise(Exposure = mean(na.omit(exposure)))


cartog_compl_expo <- left_join(cartog_prod, df_expo_5deg)

#pdf("Fig/Chap5_Fig2_Completeness_cartogram.pdf")
ggplot()+
  geom_sf(data = grid, fill = "grey60", color = NA)+
  geom_sf(data = cartog_compl_expo, aes(fill=Exposure), color = NA) +
  scale_fill_viridis_c(option="viridis")+
  theme_classic()
#dev.off()


# add relationship between exposure and completeness
expo_norm <- list()
for(i in 1:length(dl_expo_norm)){
  expo_norm[[i]] <- dl_expo_norm[[i]] %>%
    mutate(class = names(dl_expo_norm)[i]) %>%
    select(grid_id, class, contains("expo"), SR_tot_ias_log)
}
df_expo_norm <- bind_rows(expo_norm)

rel_expo_comp <- left_join(df_expo_norm, df_compl %>% st_drop_geometry()) %>%
  filter(class!="all_groups") %>%
  mutate(Class=case_when(
    class=="bird" ~ "Aves",
    class=="mam" ~ "Mammalia",
    class=="rept" ~ "Reptilia"
  )) %>%
  mutate(expo_norm = expo_max_min/max(expo_max_min))

#pdf("Fig/Chap5_Fig2b_Expo_Complet.pdf", 3, 4)
ggplot(rel_expo_comp, 
       aes(x = comp_se_product, y = expo_norm))+
  geom_point(alpha = .3, color = "grey60", size=2)+
  geom_smooth(aes(color=Class), method = "lm")+
  xlab("Completeness")+
  ylab("Exposure")+
  theme_classic()+
  theme(legend.position = "top")
#dev.off()

mod = lm(expo_norm~comp_se_product+class, data = rel_expo_comp)
summary(mod)


mod = lm(SR_tot_ias_log~comp_se_product+class, data = rel_expo_comp)
summary(mod)


plot(mod)
cor.test(rel_expo_comp$SR_tot_ias_log, rel_expo_comp$comp_se_product)
cor.test(rel_expo_comp$expo_norm, rel_expo_comp$comp_se_product)


#### Sensitivity ####


# range of completeness values
# for each taxonomic group
lapply(df_compl_sens, summary)

# set taxonomic group for exploring components of completeness
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

# final completeness = same as ed
ggplot(grid_compl_sens)+
  geom_sf(aes(fill=comp_prod), color=NA) +
  scale_fill_viridis_c(na.value = "grey80")+
  theme_classic()+
  ggtitle("Completeness")



#### Cartogram of sensitivity with completeness

# cartogram function for each group
# class="rept"

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
# pdf("Fig/Fig3_Sensit_compl_map_BMR.pdf", 9,8)
ggarrange(cb, cm, cr, legend = "top", ncol = 2, nrow=2)
# dev.off()



# Sensitivity ~ completeness 


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

sens_compl("bird")
sens_compl("mam")
sens_compl("rept")

