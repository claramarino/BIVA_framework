# Supplementary analyses: 
# 1. Relationships between normalization methods
#   - for exposure
#   - for sensitivity
#   - for VU see BIVA_Final_Vu_all_norm_methods.R
# 2. Relationships between exposure metrics 
# 3. Product vs sum of exposure components
# 4. Relationship between sensitivity and species richness

#### Load data ####

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
dl_expo_norm <- readRDS(paste0("Output/Exposure/RISK_23_expo_norm_", res, "_km"))
dl_sensit <- readRDS(paste0("Output/Sensitivity/RISK_33_sensit_norm_", res, "_km"))


#### 1 Relationships norm method ####

### For exposure ###

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

#pdf("Fig/Suppl_Fig3A_Norm_exposure.pdf", 6, 5)
ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)
#dev.off()


### For sensitivity ###

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


#pdf("Fig/Suppl_Fig3B_Norm_sensitivity.pdf", 6, 5)
ggarrange(log_max_min, log_rank,rank_max_min,
          ncol = 2, nrow = 2, common.legend = T)
#dev.off()

lapply(dl_sensit, function(x){
  round(cor(x %>% select(contains("SR")), method = "pearson"), 3)
})

#### 2 Relationships between Exposure components ####

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

#pdf("Fig/Chap5_Fig0_Rel_exposure_components.pdf", 6, 5)
ggarrange(sr_range, imp_range, imp_sr,
          ncol = 2, nrow = 2, common.legend = T)
#dev.off()


lapply(dl_metrics, function(x){
  cor(x %>% select(-grid_id), 
      method = "pearson")
})

#### Map the 3 exposure metrics to see the spatial differences

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

#pdf("Fig/Chap5_FigS3_spatial_expo_metrics.pdf", 8, 12)
ggarrange(all_sr, all_range, all_impact, ncol = 1, nrow=3,
          common.legend = T, legend = "top")
#dev.off()

#### 3 Sum vs product of exposure components #####

sum_prod <- readRDS(paste0("Output/Exposure/BIVA_10_expo_norm_", res, "_km"))


df_sum_prod <- sum_prod[["all_groups"]]


corr = cor(df_sum_prod %>% select(contains("max_min")))
p.mat <- ggcorrplot::cor_pmat(df_sum_prod %>% select(contains("max_min")))
ggcorrplot::ggcorrplot(corr, p.mat = p.mat, hc.order = TRUE,
                       type = "lower", insig = "blank", lab = T)

ggplot(df_sum_prod)+
  geom_point(aes(x=expo_max_min, y = expo_prod_max_min))+
  geom_smooth(aes(x=expo_max_min, y = expo_prod_max_min))

summary(df_sum_prod %>% select(contains("max_min")))


grid_all <- left_join(grid, df_sum_prod)

all_sum <- ggplot(grid_all)+
  geom_sf(aes(fill = expo_max_min), color = NA)+
  scale_fill_viridis_c(option="viridis", direction = -1, na.value = "grey80" )+
  theme_classic()+labs(title="Sum of components")
all_prod <- ggplot(grid_all)+
  geom_sf(aes(fill = expo_prod_max_min), color = NA)+
  scale_fill_viridis_c(option="viridis", direction = -1, na.value = "grey80" )+
  theme_classic()+labs(title="Product of components")

png("Fig/Suppl_Fig4_exposure_product_vs_sum.png", 
    height = 5, width = 9, units = "in", res=300)
ggarrange(all_sum, all_prod, ncol=2, nrow=1, legend = "top")
dev.off()

#### 4 Sensitivity and species richness #####

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

# plot final output: residuals of the linear relationship
# between sensitivity (number of native species affected by an IAS)
# and the total species richness

pdf(paste0("Fig/RISK_33_resid_iasA_SRtot_BMR_", res, ".pdf"), 8, 5)
ggarrange(b, m, r, ncol = 2, nrow=2, common.legend = T)
dev.off()
