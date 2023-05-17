# Compute sensitivity to IAS

rm(list=ls())

library(tidyverse)
library(sf)

# choose resolution 
res = "110" # 55 or 110 km
# grid
grid_terr <- readRDS(paste0("Output/RISK_32_grid_", res, "km"))

#### Load data ####

# sp presence per cell
sp_fold = "Output/Sensitivity/Cells_nat_all_bmr/"
sp_b1 <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_BIRD_110_55_chunk1"))
sp_b2 <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_BIRD_110_55_chunk2"))
sp_b3 <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_BIRD_110_55_chunk3"))
sp_b4 <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_BIRD_110_55_chunk4"))
sp_m <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_MAM_110_55"))
sp_r <- readRDS(paste0(sp_fold, "RISK_32_cells_nat_all_REPT_110_55"))

# load species id
path_poly <- "Output/Sensitivity/Polygons_native_all_bmr/"
mam_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_MAM"))
mam_IDs <- mam_nat %>% st_drop_geometry() %>%
  dplyr::mutate(objectid = 1:nrow(mam_nat)) %>%
  select(objectid, binomial)
rm(mam_nat)
b1_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD1"))
b1_IDs <- b1_nat %>% st_drop_geometry() %>%
  dplyr::mutate(objectid = 1:nrow(b1_nat)) %>%
  select(objectid, binomial)
rm(b1_nat)
b2_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD2"))
b2_IDs <- b2_nat %>% st_drop_geometry()  %>%
  dplyr::mutate(objectid = 1:nrow(b2_nat)) %>%
  select(objectid, binomial)
rm(b2_nat)
b3_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD3"))
b3_IDs <- b3_nat %>% st_drop_geometry()  %>%
  dplyr::mutate(objectid = 1:nrow(b3_nat)) %>%
  select(objectid, binomial)
rm(b3_nat)
b4_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_BIRD4"))
b4_IDs <- b4_nat %>% st_drop_geometry()  %>%
  dplyr::mutate(objectid = 1:nrow(b4_nat)) %>%
  select(objectid, binomial)
rm(b4_nat)
bird_IDs <- bind_rows(b1_IDs, b2_IDs, b3_IDs, b4_IDs)
rept_nat <- readRDS(paste0(path_poly, "RISK_32_Polygon_nat_all_REPT"))
rept_IDs <- rept_nat %>% st_drop_geometry()  %>%
  dplyr::mutate(objectid = 1:nrow(rept_nat)) %>%
  select(objectid, binomial)
rm(rept_nat)


# count total number of species
length(unique(mam_IDs$binomial))
length(unique(rept_IDs$binomial))
length(unique(bird_IDs$binomial))


# load IAS-A info

# load ias list
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")
# Load occurrence ias
in_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(in_fold, "RISK_21_grid_cells_310_IAS_", res, "km"))
# load associated natives
ias_x_native <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok")
# we keep only natives associated to the modeled ias 
native_ok <- ias_x_native %>%
  filter(ias_gbif_key %in% unique(all_sp_cells$new_key))%>%
  rename(new_key = ias_gbif_key) %>%
  # do we keep ias associated to extinct native species only?
  filter(!redlistCategory %in% c("Extinct","Extinct in the Wild")) %>%
  dplyr::select(scientificName, className, redlistCategory) %>%
  distinct() %>%
  mutate(IAS_T = if_else(
    redlistCategory %in% c("Critically Endangered","Endangered", "Vulnerable"), "IAS-T", "IAS-NT"))

length(unique(native_ok$scientificName))
table(native_ok$redlistCategory, native_ok$className)
table(native_ok$IAS_T, native_ok$className)
table(native_ok$className)
table(native_ok$IAS_T)

#### calculate species richness per cell for different groups ####

# extract cells ids for each sp 
sp_ids = b1_IDs
sp_cells = sp_b1

get_cell_ids <- function(sp_ids, sp_cells){
  df <- data.frame()
  for(i in sp_ids$objectid){
    cells = sp_cells[[i]][[paste0("cells",res)]]
    if(is_empty(cells)){
      sp_df = data.frame(
        binomial = sp_ids$binomial[i],
        cell_id =  NA)
    } else {
      sp_df = data.frame(
        binomial = sp_ids$binomial[i],
        cell_id =  cells)
    }
    df <- bind_rows(df, sp_df)
  }
  return(df %>% distinct())
}

cell_ids_m <- get_cell_ids(sp_ids = mam_IDs, sp_cells = sp_m)
saveRDS(cell_ids_m, paste0(sp_fold, "RISK_33_cell_IDs_MAM_all_nat_", res))
sr_tot_m <- cell_ids_m %>% 
  group_by(cell_id) %>%
  summarise(SR_tot_mam = n())
hist(sr_tot_m$SR_tot_mam)

cell_ids_b1 <- get_cell_ids(sp_ids = b1_IDs, sp_cells = sp_b1)
saveRDS(cell_ids_b1, paste0(sp_fold, "RISK_33_cell_IDs_BIRD1_all_nat_", res))
cell_ids_b2 <- get_cell_ids(sp_ids = b2_IDs, sp_cells = sp_b2)
saveRDS(cell_ids_b2, paste0(sp_fold, "RISK_33_cell_IDs_BIRD2_all_nat_", res))
cell_ids_b3 <- get_cell_ids(sp_ids = b3_IDs, sp_cells = sp_b3)
saveRDS(cell_ids_b3, paste0(sp_fold, "RISK_33_cell_IDs_BIRD3_all_nat_", res))
cell_ids_b4 <- get_cell_ids(sp_ids = b4_IDs, sp_cells = sp_b4)
saveRDS(cell_ids_b4, paste0(sp_fold, "RISK_33_cell_IDs_BIRD4_all_nat_", res))
cell_ids_r <- get_cell_ids(sp_ids = rept_IDs, sp_cells = sp_r)
saveRDS(cell_ids_r, paste0(sp_fold, "RISK_33_cell_IDs_REPT_all_nat_", res))


cell_ids <- list(
  mam = cell_ids_m,
  bird = bind_rows(cell_ids_b1, cell_ids_b2, cell_ids_b3, cell_ids_b4) %>% distinct(),
  rept = cell_ids_r
)

sr_tot <- lapply(cell_ids, function(x){
  x %>% 
    group_by(cell_id) %>%
    summarise(SR_tot = n())
})

saveRDS(sr_tot, paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_", res))

sr_ias_a <- lapply(cell_ids, function(x){
  x %>% 
    filter(binomial %in% native_ok$scientificName) %>%
    group_by(cell_id) %>%
    summarise(SR_ias_a = n())
})

sr_ias_t <- lapply(cell_ids, function(x){
  x %>% 
    filter(binomial %in% native_ok$scientificName[native_ok$IAS_T=="IAS-T"]) %>%
    group_by(cell_id) %>%
    summarise(SR_ias_t = n())
})

sr_ias_nt <- lapply(cell_ids, function(x){
  x %>% 
    filter(binomial %in% native_ok$scientificName[native_ok$IAS_T=="IAS-NT"]) %>%
    group_by(cell_id) %>%
    summarise(SR_ias_nt = n())
})

sr_all <- mapply(left_join, 
              mapply(full_join, sr_tot, sr_ias_a, SIMPLIFY = FALSE), 
              mapply(full_join, sr_ias_t, sr_ias_nt, SIMPLIFY = FALSE), 
              SIMPLIFY = FALSE)
colnames(sr_all$mam)

sr_all <- lapply(sr_all, function(x){
  x %>% replace_na(list(SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0)) %>%
    filter(!is.na(cell_id))
})

saveRDS(sr_all, paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_ias_a_t_nt_", res))


##### Calculate sensitivity metrics with normalization #####

sr_all <- readRDS(paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_ias_a_t_nt_", res))

grid_terr$cell_id = 1:nrow(grid_terr)


dl_all_norm <- lapply(sr_all, function(x){
  
  #x=sr_all$mam # for test
  
  #------ initialisation
  dl_max_min  <- dl_log <- dl_rank <- x
  
  # max min linear
  maxcol <- apply(x, 2, max)
  mincol <- apply(x, 2, min)
  # log transformed
  maxlog <- apply(x, 2, function(x){max(log(x+1))})
  minlog <- apply(x, 2, function(x){min(log(x+1))})
  # ranks
  x_rank <- x %>%
    dplyr::mutate_at(vars(-cell_id), dense_rank)
  maxrank <- apply(x_rank, 2, max)
  minrank <- apply(x_rank, 2, min)
  
  #------ normalization
  for (i in 2:length(maxcol)){
    # max min
    dl_max_min[,i] <- (x[,i]-mincol[i])/(maxcol[i]-mincol[i])
    # log-transformed
    dl_log[,i] <- ((log(x[,i]+1)-minlog[i]))/
      (maxlog[i]-minlog[i])
    # ranks
    dl_rank[,i] <- (x_rank[,i]-minrank[i])/(maxrank[i]-minrank[i])
  }
  
  #------ correct names
  colnames(dl_max_min)[2:ncol(dl_max_min)] <- paste0(colnames(x)[2:ncol(x)], "_max_min")
  colnames(dl_log)[2:ncol(dl_log)] <- paste0(colnames(x)[2:ncol(x)], "_log")
  colnames(dl_rank)[2:ncol(dl_rank)] <- paste0(colnames(x)[2:ncol(x)], "_rank")
  
  # bind all tables
  all_norm <- list(
    dl_max_min, dl_log, dl_rank) %>% 
    reduce(left_join, by="cell_id")
  
  return(all_norm)
  
})

lapply(dl_all_norm, summary)

# keep only metrics related to ias-a
dl_all_norm <- lapply(dl_all_norm, function(x){
  x %>% dplyr::select(cell_id, contains("SR_ias_a"))
})

# save sensitivity table
saveRDS(dl_all_norm, paste0("Output/Sensitivity/RISK_33_sensit_norm_", res, "_km"))




############### Tests de plots ##################

grid_mam <- left_join(grid_terr, sr_all$mam) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))


ggplot(grid_mam, aes(x=SR_tot, y = SR_ias_a))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(grid_mam, aes(x=log(SR_tot+1), y = log(SR_ias_a+1)))+
  geom_point()+
  geom_smooth(method = "lm")

lm_mam <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_mam)
summary(lm_mam)
grid_mam$res = residuals(lm_mam)
hist(grid_mam$res)

ggplot(grid_mam)+
  geom_sf(aes(fill=SR_ias_a), color=NA) +
  scale_fill_viridis_c(option = "inferno", direction = -1)+
  theme_classic()+
  ggtitle("Mammalia")

m<- ggplot(grid_mam)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Mammalia")

m

library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE) 
RColorBrewer::brewer.pal(name = "RdYlBu", n=11)


# birds

grid_bird <- left_join(grid_terr, sr_all$bird) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

ggplot(grid_bird, aes(x=SR_tot, y = SR_ias_a))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(grid_bird, aes(x=log(SR_tot+1), y = log(SR_ias_a+1)))+
  geom_point()+
  geom_smooth(method = "lm")

lm_bird <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_bird)
summary(lm_bird)
grid_bird$res = residuals(lm_bird)
hist(grid_bird$res)

ggplot(grid_bird)+
  geom_sf(aes(fill=SR_ias_a), color=NA) +
  scale_fill_gradient(
    low = "white", 
    high = "red")

b<-ggplot(grid_bird)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Aves")
b
ggplot(grid_bird)+
  geom_sf(aes(fill=SR_ias_a), color=NA)  +
  scale_fill_viridis_c(option = "inferno", direction = -1)+
  theme_classic()+
  ggtitle("Aves")


#### reptiles


grid_rept <- left_join(grid_terr, sr_all$rept) %>%
  replace_na(list(SR_tot = 0, SR_ias_a = 0, SR_ias_t = 0, SR_ias_nt = 0))

ggplot(grid_rept)+
  geom_sf(aes(fill=SR_tot), color=NA) +
  scale_fill_gradient(
    low = "white", 
    high = "red")

ggplot(grid_rept, aes(x=SR_tot, y = SR_ias_a))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(grid_rept, aes(x=log(SR_tot+1), y = log(SR_ias_a+1)))+
  geom_point()+
  geom_smooth(method = "lm")

lm_rept <- lm(log(SR_ias_a+1)~log(SR_tot+1), data = grid_rept)
summary(lm_rept)
grid_rept$res = residuals(lm_rept)
hist(grid_rept$res)

ggplot(grid_rept)+
  geom_sf(aes(fill=SR_ias_a), color=NA) +
  scale_fill_viridis_c(option = "inferno", direction = -1)+
  theme_classic()+
  ggtitle("Reptilia")

r <-  ggplot(grid_rept)+
  geom_sf(aes(fill=res), color=NA) +
  scale_fill_gradient2(low = "#313695", mid = "#FFFFBF", high = "#D73029")+
  theme_classic()+
  ggtitle("Reptilia")

r







###############

sr_b <- bind_rows(sr_b1, sr_b2, sr_b3, sr_b4) %>%
  group_by(x, y) %>%
  summarise(SR_tot = sum (SR_tot))
rm(sr_b1, sr_b2, sr_b3, sr_b4)


sr_bmr <- inner_join(
  area,
  bind_rows(sr_b, sr_m, sr_r) %>%
    group_by(x, y) %>%
    summarise(SR_tot_bmr = sum (SR_tot)))


sr_bmr_area <- sr_bmr %>%
  mutate(SR_tot_bmr_area = SR_tot_bmr/area)
ggplot(sr_bmr_area, aes(y, SR_tot_bmr_area))+
  geom_smooth()

#### Compute metrics ####


# for each metric, compute a SR total across BMR

# ias a all
sr_iasa_bmr <- bind_rows(
  sr_iasa_b %>% rename(SR_tot = SR_iasa_b), 
  sr_iasa_m %>% rename(SR_tot = SR_iasa_m), 
  sr_iasa_r %>% rename(SR_tot = SR_iasa_r)) %>%
  group_by(x, y) %>%
  summarise(SR_iasa_bmr = sum(SR_tot))
# ias a mod
sr_iasa_mod_bmr <- bind_rows(
  sr_iasa_mod_b %>% rename(SR_tot = SR_iasa_mod_b), 
  sr_iasa_mod_m %>% rename(SR_tot = SR_iasa_mod_m), 
  sr_iasa_mod_r %>% rename(SR_tot = SR_iasa_mod_r)) %>%
  group_by(x, y) %>%
  summarise(SR_iasa_mod_bmr = sum(SR_tot))

# ias t all
sr_iast_bmr <- bind_rows(
  sr_iast_b %>% rename(SR_tot = SR_iast_b), 
  sr_iast_m %>% rename(SR_tot = SR_iast_m), 
  sr_iast_r %>% rename(SR_tot = SR_iast_r)) %>%
  group_by(x, y) %>%
  summarise(SR_iast_bmr = sum(SR_tot))
# ias t mod
sr_iast_mod_bmr <- bind_rows(
  sr_iast_mod_b %>% rename(SR_tot = SR_iast_mod_b), 
  sr_iast_mod_m %>% rename(SR_tot = SR_iast_mod_m), 
  sr_iast_mod_r %>% rename(SR_tot = SR_iast_mod_r)) %>%
  group_by(x, y) %>%
  summarise(SR_iast_mod_bmr = sum(SR_tot))


# Bind all metrics in one df

dfs <- list(
  sr_bmr, sr_iasa_bmr, sr_iast_bmr, sr_iasa_mod_bmr, sr_iast_mod_bmr,
  sr_b %>% rename(SR_tot_b = SR_tot), sr_iasa_b, sr_iast_b, sr_iasa_mod_b, sr_iast_mod_b,
  sr_m %>% rename(SR_tot_m = SR_tot), sr_iasa_m, sr_iast_m, sr_iasa_mod_m, sr_iast_mod_m,
  sr_r %>% rename(SR_tot_r = SR_tot), sr_iasa_r, sr_iast_r, sr_iasa_mod_r, sr_iast_mod_r)


sensit_tbl <- plyr::join_all( dfs, by = c("x","y"), type = "full")
str(sensit_tbl)

# replace NA with 0
sensit_tbl[is.na(sensit_tbl)] <- 0
str(sensit_tbl)


# save output 
saveRDS(sensit_tbl, paste0("Output/Sensitivity/RISK_33_sensitivity_table_r", deg))



##### Explore relations between metrics #####

deg = "01"
sensit_tbl <- readRDS(paste0("Output/Sensitivity/RISK_33_sensitivity_table_r", deg))

summary(sensit_tbl)
rm(area)
# correct all metrics by area 
sensit_tbl_corr <- sensit_tbl %>%
  mutate_at(vars(contains('SR_')), list(by_area=~.*100/area))


# normalize sensitivity values between 0-1 to compare cells 
# 3 methods for rescaling variables

# max min linear

maxcol <- apply(sensit_tbl_corr, 2, max)
mincol <- apply(sensit_tbl_corr, 2, min)

sensit_tbl_norm_lin <- sensit_tbl_corr
for (i in 3:length(maxcol)){
  sensit_tbl_norm_lin[,i] <- (sensit_tbl_norm_lin[,i]-mincol[i])/(maxcol[i]-mincol[i])
}

summary(sensit_tbl_norm_lin)
hist(sensit_tbl_norm_lin$SR_iasa_bmr)
hist(sensit_tbl_norm_lin$SR_iasa_bmr_by_area)


#log transform and max min linear
maxlog <- apply(sensit_tbl_corr, 2, function(x){max(log(x+1))})
minlog <- apply(sensit_tbl_corr, 2, function(x){min(log(x+1))})
sensit_tbl_norm_log <- sensit_tbl_corr
for (i in 3:length(maxlog)){
  sensit_tbl_norm_log[,i] <- (log(sensit_tbl_norm_log[,i]+1)-minlog[i])/
    (maxlog[i]-minlog[i])
}
summary(sensit_tbl_norm_log)
hist(sensit_tbl_norm_log$SR_iasa_bmr)
hist(sensit_tbl_norm_log$SR_iasa_bmr_by_area)


# cumulative distribution
var_rank <- sensit_tbl_corr %>% 
  dplyr::select(-c(x, y, area)) %>% 
  mutate_all(dense_rank)
sensit_tbl_rank <- bind_cols(sensit_tbl_corr %>% dplyr::select(x,y), var_rank)

maxrank <- apply(sensit_tbl_rank, 2, max)
minrank <- apply(sensit_tbl_rank, 2, min)

for (i in 3:length(maxrank)){
  sensit_tbl_rank[,i] <- (sensit_tbl_rank[,i]-minrank[i])/(maxrank[i]-minrank[i])
}

summary(sensit_tbl_rank)



saveRDS(sensit_tbl_norm_log, 
        paste0("Output/Sensitivity/RISK_33_sensitivity_norm_log_r", deg))
saveRDS(sensit_tbl_norm_lin, 
        paste0("Output/Sensitivity/RISK_33_sensitivity_norm_lin_r", deg))
saveRDS(sensit_tbl_rank, 
        paste0("Output/Sensitivity/RISK_33_sensitivity_norm_rank_r", deg))





# correlations betweens normalized var
library(ggcorrplot)
mcor <- cor(sensit_tbl_norm_log)
ggcorrplot(mcor[-c(1,2), -c(1,2)])

# map normalized sensitivity

#bmr IAS-A
snorm <- ggplot(data = sensit_tbl_norm_log) +
  geom_raster(aes(x = x, y = y, fill = SR_iasa_bmr_by_area)) +
  scale_fill_gradient(
    low = "white", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
snorm



ggplot(sensit_tbl_norm_lin)+
  geom_smooth(aes(y, SR_iasa_bmr_by_area))+
  geom_smooth(aes(y, SR_iasa_bmr), color = "red")


#bmr IAS-T
snorm <- ggplot(data = sensit_tbl_norm_log) +
  geom_raster(aes(x = x, y = y, fill = SR_iast_bmr)) +
  scale_fill_gradient(
    low = "white", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
snorm


ggplot(sensit_tbl_norm) +
  geom_point(aes(x=SR_iast_bmr, y=SR_iasa_bmr))

ggplot(sensit_tbl_norm) +
  geom_point(aes(x=SR_iasa_m, y=SR_iasa_r))
ggplot(sensit_tbl_norm) +
  geom_point(aes(x=SR_iasa_b, y=SR_iasa_r))
ggplot(sensit_tbl_norm) +
  geom_point(aes(x=SR_iasa_m, y=SR_iasa_b))

hist(sensit_tbl_norm$SR_iasa_bmr)
hist(sensit_tbl_norm$SR_iasa_b)
hist(sensit_tbl_norm$SR_iasa_m)
hist(sensit_tbl_norm$SR_iasa_r)












# try correlations between raw variables 
mcor <- cor(sensit_tbl)
ggcorrplot(mcor[-c(1,2), -c(1,2)])

norm_sens <- sensit_tbl %>% 
  mutate(ias_a_norm = SR_iasa_m/SR_tot_m,
         ias_t_norm = SR_iast_m/SR_tot_m)

summary(norm_sens)

SR_bmr_a <- ggplot(data = norm_sens %>% filter(ias_a_norm>0)) +
  geom_raster(aes(x = x, y = y, fill = ias_a_norm)) +
  scale_fill_gradient(
    low = "grey", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_bmr_a
SR_bmr_a <- ggplot(data = norm_sens %>% filter(ias_t_norm>0)) +
  geom_raster(aes(x = x, y = y, fill = ias_t_norm)) +
  scale_fill_gradient(
    low = "grey", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_bmr_a

SR_bmr_a <- ggplot(data = sensit_tbl %>% filter(SR_iasa_bmr>0)) +
  geom_raster(aes(x = x, y = y, fill = SR_iasa_bmr)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_bmr_a

SR_bmr_t <- ggplot(data = sensit_tbl %>% filter(SR_iast_bmr>0)) +
  geom_raster(aes(x = x, y = y, fill = SR_iast_bmr)) +
  scale_fill_gradient(
    low = "green", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic()
SR_bmr_t


