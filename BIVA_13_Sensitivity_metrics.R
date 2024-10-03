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
saveRDS(dl_all_norm, paste0("Output/Sensitivity/BIVA_13_sensit_norm_", res, "_km"))

openxlsx::write.xlsx(dl_all_norm, file = paste0("Output/Sensitivity/BIVA_13_sensit_norm_", res, "_km.xlsx"))
sf::st_write(grid_terr, paste0("Output/Sensitivity/Grid_", res, "_km.shp"))

