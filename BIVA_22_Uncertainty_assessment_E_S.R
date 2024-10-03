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

# check completeness per geographic zones
# for each continent

# import a world countries map:
countries <- geodata::world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
# import a table with country codes and continents:
cntry_codes <- geodata::country_codes()

# project to good CRS
# make the geometries valid for intersection bewteen grid and countries
countries <- sf::st_as_sf(countries)
df_compl <- sf::st_transform(df_compl, sf::st_crs(countries))
sf::st_crs(countries) == sf::st_crs(df_compl)
sum(sf::st_is_valid(countries))
countries <- sf::st_make_valid(countries)
df_compl <- sf::st_make_valid(df_compl)

sum(sf::st_is_valid(countries))
sum(sf::st_is_valid(df_compl))

compl_count <- sf::st_intersects(countries, df_compl)
sum(unlist(lapply(compl_count, is_empty)))
 for(i in 1:length(compl_count)){
   compl_count[[i]] <- data.frame(
     cell_id = compl_count[[i]],
     ISO3 = countries$GID_0[i]
   )
 }
compl_df <- bind_rows(compl_count)

compl_df_all <- left_join(left_join(compl_df, cntry_codes),
                          left_join(grid %>% sf::st_drop_geometry(), df_compl %>% sf::st_drop_geometry()))

unique(compl_df_all$UNREGION1)
unique(compl_df_all$UNREGION2)
unique(compl_df_all$continent)

mean_compl <- compl_df_all %>%
  group_by(continent) %>%
  summarize(mean_compl = mean(comp_se_product, na.rm = T),
            mean_se = mean(norm_dens, na.rm = T),
            mean_ed = mean(DB_and_GRIIS, na.rm = T))

png("Fig/Suppl_Fig5_Completeness_expo_Continents.png",
    height = 5 , width = 5, units = "in", res=300)
ggplot(compl_df_all)+
  geom_point(aes(x=comp_se_product, y = continent, color = continent), 
             alpha = .2, position = "jitter")+
  geom_boxplot(aes(x=comp_se_product, y = continent), fill=NA, outlier.shape = NA)+
  geom_point(data = mean_compl, aes(x=mean_compl, y = continent), shape = 23)+
  xlab("Completeness of exposure to IAS")+ ylab("Continent")+
  theme_bw()
dev.off()


png("Fig/Suppl_Fig5_Completeness_expo_Regions.png",
    height = 7 , width = 6, units = "in", res=300)
ggplot(compl_df_all)+
  geom_point(aes(x=comp_se_product, y = UNREGION1, color = continent), 
             alpha = .4, position = "jitter")+
  geom_boxplot(aes(x=comp_se_product, y = UNREGION1), fill=NA, outlier.shape = NA)+
  xlab("Completeness of exposure to IAS")+ ylab("UN Region")+
  theme_bw()
dev.off()

mean_compl_reg <- compl_df_all %>%
  group_by(UNREGION1) %>%
  summarize(mean_compl = mean(comp_se_product, na.rm = T),
            med_compl = median(comp_se_product, na.rm = T),
            mean_se = mean(norm_dens, na.rm = T),
            mean_ed = mean(DB_and_GRIIS, na.rm = T))

ggplot(compl_df_all)+
  geom_point(aes(x=norm_dens, y = continent, color = continent), 
             alpha = .2, position = "jitter")+
  geom_boxplot(aes(x=norm_dens, y = continent), fill=NA, outlier.shape = NA)+
  geom_point(data = mean_compl, aes(x=mean_se, y = continent), shape = 23)+
  xlab("Sampling effort")+ ylab("Continent")+
  theme_bw()

ggplot(compl_df_all)+
  geom_point(aes(x=DB_and_GRIIS, y = continent, color = continent), 
             alpha = .2, position = "jitter")+
  geom_boxplot(aes(x=DB_and_GRIIS, y = continent), fill=NA, outlier.shape = NA)+
  geom_point(data = mean_compl, aes(x=mean_ed, y = continent), shape = 23)+
  xlab("Effective detection")+ ylab("Continent")+
  theme_bw()

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


# check completeness per geographic zones
# for mammals and reptiles because birds are well informed

#select the group
group = "rept" # mam or rept

# load the intersection computed for exposure (L133 to 163)

compl_df_sens <- left_join(compl_df_all, df_compl_sens[[group]])

mean_compl <- compl_df_sens %>%
  group_by(continent) %>%
  summarize(mean_compl = mean(comp_prod, na.rm = T))

png(paste0("Fig/Suppl_Fig5_Completeness_sensi_", group, "_Continents.png"),
    height = 5 , width = 5, units = "in", res=300)
ggplot(compl_df_sens)+
  geom_point(aes(x=comp_prod, y = continent, color = continent), 
             alpha = .2, position = "jitter")+
  geom_boxplot(aes(x=comp_prod, y = continent), fill=NA, outlier.shape = NA)+
  geom_point(data = mean_compl, aes(x=mean_compl, y = continent), shape = 23)+
  xlab(paste0("Completeness of exposure to IAS (", group, ")"))+ ylab("Continent")+
  theme_bw()
dev.off()


png(paste0("Fig/Suppl_Fig5_Completeness_sensi_", group, "_Regions.png"),
    height = 7 , width = 6, units = "in", res=300)
ggplot(compl_df_sens)+
  geom_point(aes(x=comp_prod, y = UNREGION1, color = continent), 
             alpha = .4, position = "jitter")+
  geom_boxplot(aes(x=comp_prod, y = UNREGION1), fill=NA, outlier.shape = NA)+
  xlab(paste0("Completeness of exposure to IAS (", group, ")"))+ ylab("UN Region")+
  theme_bw()
dev.off()

####### Include a cutoff on completeness for VU ######

########## Tests 


hist(df_compl$comp_se_product)
quantile(df_compl$comp_se_product)

expo_comp <- left_join(dl_expo_norm$all_groups %>% dplyr::select(grid_id, expo_max_min), 
                       df_compl) %>% sf::st_as_sf()

hist(expo_comp$comp_se_product)

# adopt a cutoff based on values' distribution
quantile(expo_comp$comp_se_product, probs = c(1/3, .5, 2/3, 1))

expo_comp_cut <- expo_comp %>%
  mutate(comp.3 = if_else(comp_se_product>0.23051245 , 1, 0), # 1/3 of values with the best completeness
         comp.5 = if_else(comp_se_product>0.1319, 1, 0))%>% # 50% of values with the best completeness
  mutate(Exposure = expo_max_min/max(expo_max_min),
         comp01 = comp_se_product/max(comp_se_product)) %>%
  mutate(comp01.6 = if_else(comp01>= 2/3 , 1, 0),
         comp01.3 = if_else(comp01>= 1/3 , 1, 0))

class(expo_comp_cut)

ggplot(expo_comp_cut)+
  geom_sf(aes(fill = Exposure), color = NA) +
  scale_fill_viridis_c(option="inferno", direction = -1)+
  theme_classic()

ggplot(expo_comp_cut)+
  geom_sf(aes(fill = Exposure, alpha = comp01), color = NA) +
  scale_fill_viridis_c(option="inferno", direction = -1)+
  theme_classic()

ggplot(expo_comp_cut)+
  geom_sf(aes(fill = Exposure, alpha = comp.5), color = NA) +
  scale_fill_viridis_c(option="inferno", direction = -1)+
  theme_classic()
ggplot(expo_comp_cut)+
  geom_sf(aes(fill = Exposure, alpha = comp.3), color = NA) +
  scale_fill_viridis_c(option="inferno", direction = -1)+
  theme_classic()

hist(expo_comp_cut$comp01)

ggplot(expo_comp_cut)+
  geom_sf(aes(fill = Exposure, alpha = comp01.6), color = NA) +
  scale_fill_viridis_c(option="inferno", direction = -1)+
  theme_classic()

class(expo_comp)
class(df_compl)


View(df_compl_sens$mam)

######################## Clean output


ntot = nrow(grid)

# percentage of cells with at least one of the 304 IAS
lapply(dl_expo_norm, function(x){
  nrow(x)/ntot
})

# percentage of cells with at least one sensitive native species to the 304 IAs
lapply(df_compl_sens, function(x){
  nrow(x %>% filter(SR_ias_a>0))/17257
})

class = "bird"
cut = 0.25 # 0.25 = 0.5*0.5; 0.5625 = 0.75*0.75


#### Exposure

plot_expo_cut <- function(class, cut){
  # class is bird, mam, rept
  if(!class %in% c("bird","mam","rept")){
    stop("class must be 'bird', 'mam', or 'rept")}
  
  # cut is numeric: .3, .5, .6
  
  dat <- dl_expo_norm[[class]] 
  dat_sf <- left_join(grid, dat) %>%
    select(grid_id, cell_id, contains("expo")) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min))
  sf_cut <- left_join(dat_sf, df_compl %>% sf::st_drop_geometry()) %>%
    #mutate(comp01 = comp_se_product/max(comp_se_product)) %>%
    mutate(comp_cut = if_else(comp_se_product>= cut , 1, 0))
  p <- ggplot(sf_cut)+
    geom_sf(data=grid, fill="grey90", color = NA)+
    geom_sf(aes(fill = Exposure, alpha = comp_cut), color = NA) +
    scale_fill_viridis_c(option="inferno", direction = -1)+
    scale_alpha_continuous(range=c(0,1), limits=c(0.5,1), na.value = 0)+
    theme_classic()
  return(p)
}


eb <- plot_expo_cut("bird", 0)
eb1 <- plot_expo_cut("bird", 1/9)
eb2 <- plot_expo_cut("bird", 4/9)

em <- plot_expo_cut("mam", 0)
em1 <- plot_expo_cut("mam", 1/9)
em2 <- plot_expo_cut("mam", 4/9)

er <- plot_expo_cut("rept", 0)
er1 <- plot_expo_cut("rept", 1/9)
er2 <- plot_expo_cut("rept", 4/9)

e <- ggarrange(eb, em, er,
               ncol = 1, nrow = 3, legend = "none")
e1 <- ggarrange(eb1, em1, er1,
               ncol = 1, nrow = 3, legend = "none")
e2 <- ggarrange(eb2, em2, er2,
               ncol = 1, nrow = 3, legend = "none")


# Save final figure
pdf("Fig/SupplFig4_Exposure_with_cutoffs.pdf", 8,6)
ggarrange(e, e1, e2, ncol = 3, common.legend = F)
dev.off()
png("Fig/SupplFig4_Exposure_with_cutoffs.png", width = 1600, height = 1200, res = 300 )
ggarrange(e, e1, e2, ncol = 3, common.legend = F)
dev.off()


legend_e <- cowplot::get_legend(er1+ theme(legend.direction = "horizontal"))


#### Sensitivity

plot_sensib_cut <- function(class, cut){
  dat <- dl_sensit[[class]] 
  dat_sf <- left_join(grid, dat) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  dat_compl <- left_join(grid %>% sf::st_drop_geometry(), df_compl_sens[[class]]) 
  sf_cut <- left_join(dat_sf, dat_compl) %>%
    mutate(comp_cut = if_else(comp_prod>= cut , 1, 0))
  p <- ggplot(sf_cut)+
    geom_sf(data=grid, fill="grey90", color = NA)+
    geom_sf(aes(fill = Sensitivity, alpha = comp_cut), color = NA) +
    scale_fill_viridis_c(option="mako", direction = -1)+
    scale_alpha_continuous(range=c(0,1), limits=c(0.5,1), na.value = 0)+
    theme_classic()
  return(p)
}

sb <- plot_sensib_cut("bird", 0)
sb1 <- plot_sensib_cut("bird", 1/9)
sb2 <- plot_sensib_cut("bird", 4/9)

sm <- plot_sensib_cut("mam", 0)
sm1 <- plot_sensib_cut("mam", 1/9)
sm2 <- plot_sensib_cut("mam", 4/9)

sr <- plot_sensib_cut("rept", 0)
sr1 <- plot_sensib_cut("rept", 1/9)
sr2 <- plot_sensib_cut("rept", 4/9)

s <- ggarrange(sb, sm, sr,
               ncol = 1, nrow = 3, legend = "none")
s1 <- ggarrange(sb1, sm1, sr1,
               ncol = 1, nrow = 3, legend = "none")
s2 <- ggarrange(sb2, sm2, sr2,
               ncol = 1, nrow = 3, legend = "none")


# Save final figure
pdf("Fig/SupplFig4_Sensitivity_with_cutoffs.pdf", 8,6)
ggarrange(s, s1, s2, ncol = 3, common.legend = F)
dev.off()
png("Fig/SupplFig4_Sensitivity_with_cutoffs.png", width = 1600, height = 1200, res = 300 )
ggarrange(s, s1, s2, ncol = 3, common.legend = F)
dev.off()

legend_s <- cowplot::get_legend(sr + theme(legend.direction = "horizontal"))


# Define low to high exposure and sensitivity
# merge the two information
# get vulnerability categories: VL, L, M, H, VH
plot_vu_cut <- function(class, cut){
  dat1 <- left_join(dl_sensit[[class]], df_compl_sens[[class]]) %>%
    mutate(comp_cut_s = if_else(comp_prod>= cut , 1, 0)) 
  dat2 <- left_join(dl_expo_norm[[class]], df_compl %>% sf::st_drop_geometry()) %>%
    mutate(comp_cut_e = if_else(comp_se_product>= cut , 1, 0))
  
  dat_sf <- left_join(left_join(grid, dat1), dat2) %>%
    replace_na(list(expo_max_min = 0, expo_log = 0, expo_rank = 0)) %>%
    mutate(Exposure = expo_max_min/max(expo_max_min)) %>%
    replace_na(list(SR_ias_a_max_min = 0, SR_ias_a_log = 0, SR_ias_a_rank = 0))%>%
    rename(Sensitivity = SR_ias_a_max_min)
  
  data <- biscale::bi_class(dat_sf, x = Exposure, y = Sensitivity, 
                   style = "fisher", dim = 3)
  
  data_vu <- data %>%
    mutate(vu = case_when(
      bi_class=="1-1" ~ "VL",
      bi_class %in% c("1-2", "2-1") ~ "L",
      bi_class %in% c("1-3", "3-1", "2-2") ~ "M",
      bi_class %in% c("3-2", "2-3") ~ "H",
      bi_class=="3-3" ~ "VH")) %>%
    mutate(vu = factor(vu, ordered = T, levels = c("VL","L","M","H","VH"))) %>%
    mutate(vu_cut=comp_cut_s+comp_cut_e)
  
  table(data_vu$vu)
  ggplot(data_vu)+
    geom_sf(data=grid, fill="grey90", color = NA)+
    geom_sf(aes(fill=vu, alpha = vu_cut), color=NA) +
    scale_fill_manual(values = c("VL" = "#FEECEB", 
                                 "L" = "#FCC7C3", 
                                 "M" = "#F88279", 
                                 "H" = "#D63B2F", 
                                 "VH" = "#7A221B"))+
    scale_alpha_continuous(range=c(0,1), limits=c(1.5,2), na.value = 0)+ #limits=c(1.5,2) if threshold harder
    theme_classic()
}



vb <- plot_vu_cut("bird", 0)
vb1 <- plot_vu_cut("bird", 1/9)
vb2 <- plot_vu_cut("bird", 4/9)

vm <- plot_vu_cut("mam", 0)
vm1 <- plot_vu_cut("mam", 1/9)
vm2 <- plot_vu_cut("mam", 4/9)

vr <- plot_vu_cut("rept", 0)
vr1 <- plot_vu_cut("rept", 1/9)
vr2 <- plot_vu_cut("rept", 4/9)


ggarrange(vb, vb1, vb2, 
                ncol= 3, nrow=1, legend = "none")

# vulnerability figure with different cutoffs
vu <- ggarrange(vb, vm, vr, 
               ncol= 1, nrow=3, legend = "none")
vu1 <- ggarrange(vb1, vm1, vr1, 
                ncol= 1, nrow=3, legend = "none")
vu2 <- ggarrange(vb2, vm2, vr2, 
                 ncol= 1, nrow=3, legend = "none")

# Save final figure
pdf("Fig/Fig4_Vu_with_cutoffs.pdf", 8,6)
ggarrange(vu, vu1, vu2, ncol = 3, common.legend = F)
dev.off()


#save legend for VU color code
legend_vu <- cowplot::get_legend(vr + theme(legend.direction = "horizontal"))
pdf("Fig/Fig2_Legend_vu.pdf")
grid.newpage()
grid.draw(legend_vu)
dev.off()

