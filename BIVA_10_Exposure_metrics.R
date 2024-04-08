# compute & normalize exposure metrics
# save table for fig & analysis with sensitivity
# separate by victim group

rm(list=ls())

library(tidyverse)


###### Data loading #####

# select resolution
res = "110" # 55km or 110km

# Load occurrence ias
in_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(in_fold, "RISK_21_grid_cells_310_IAS_", res, "km"))
# cb de cellules avec au moins 1 IAS ?
length(unique(all_sp_cells$grid_id))
range(all_sp_cells$grid_id)

# load ias list
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")

# load ias_x_native table
ias_x_native <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok")

# do all ias in cells have associated natives in list ?
sum(unique(all_sp_cells$new_key) %in% unique(ias_x_native$ias_gbif_key))
# yes (logical but cool)


###### Native data cleaning #####

# we keep only natives associated to the modeled ias 
ias_x_native_ok <- ias_x_native %>%
  filter(ias_gbif_key %in% unique(all_sp_cells$new_key))%>%
  rename(new_key = ias_gbif_key) %>%
  # do we keep ias associated to extinct native species only?
  filter(!redlistCategory %in% c("Extinct","Extinct in the Wild")) %>%
  dplyr::select(new_key, scientificName, className, redlistCategory) %>%
  distinct()

sum(unique(all_sp_cells$new_key) %in% unique(ias_x_native_ok$new_key))


###### Compute cell related metrics #####

# combine for calculating cell related metrics
df_all <- left_join(all_sp_cells,
                    ias_x_native_ok,
                    relationship = "many-to-many")

str(df_all)
unique(df_all$className)

unique(df_all %>% filter(is.na(className)) %>% pull(new_species))
# create a list for computing metrics for all tetrapods
# and for each class separately
dl <- list(
  all_groups = df_all,
  bird = df_all %>% filter(className=="AVES"),
  mam = df_all %>% filter(className=="MAMMALIA"),
  rept = df_all %>% filter(className=="REPTILIA")
)

nrow(dl$bird)+nrow(dl$mam)+nrow(dl$rept) == nrow(dl$all_groups)

# create function for computing each metric

# function for counting ias species richness per cell
calc_sr_ias <- function(x){
  # x=dl$mam # for testing function
  out <- x %>%
    dplyr::select(grid_id, new_key) %>%
    distinct() %>%
    group_by(grid_id) %>%
    summarise(SR_tot_ias = n())
  return(out)
}

# function for ias range

calc_range <- function(x){
  # x=dl$bird
  sp_range <- x %>% 
    dplyr::select(grid_id, new_key) %>%
    distinct() %>%
    group_by(new_key) %>%
    summarise(range = n())
  out <- left_join(x %>% dplyr::select(grid_id, new_key),
                   sp_range) %>%
    distinct()%>%
    group_by(grid_id) %>%
    summarise(med_range = median(range))
return(out)
}

# function for impact breadth 
# median number of ias-a per ias per cell
calc_ib <- function(x){
  # x=dl$bird
  nb_ias_a <- x %>% 
    dplyr::select(grid_id, new_key, scientificName) %>%
    distinct() %>%
    group_by(grid_id, new_key) %>%
    summarise(ias_a = n()) #or log(n) ???
  out <- nb_ias_a %>%
    distinct()%>%
    group_by(grid_id) %>%
    summarise(med_ib = median(ias_a))
  return(out)
}


# apply functions to the list with all taxa and 1 df for each victim group

dl_metrics <- lapply(dl, function (df){
  met_to_join <- list(
    calc_sr_ias(df),
    calc_range(df),
    calc_ib(df))
  out <- met_to_join %>% reduce(left_join, by="grid_id")
  return(out)
})

# check correlations between raw metrics
lapply(dl_metrics, function(x){cor(x %>% dplyr::select(-grid_id))})



##### Normalize exposure metrics #####



dl_all_norm <- lapply(dl_metrics, function(x){
  
  #x=dl_metrics$rept # for test
  
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
    dplyr::mutate_at(vars(-grid_id), dense_rank)
  maxrank <- apply(x_rank, 2, max)
  minrank <- apply(x_rank, 2, min)
  
  #------ normalization
  for (i in 2:length(maxcol)){
    # max min
    dl_max_min[,i] <- (x[,i]-mincol[i]+0.001)/(maxcol[i]-mincol[i]+0.001)
    # log-transformed
    dl_log[,i] <- ((log(x[,i]+1)-minlog[i])+0.001)/
      (maxlog[i]-minlog[i] + 0.001)
    # ranks
    dl_rank[,i] <- (x_rank[,i]-minrank[i]+0.001)/(maxrank[i]-minrank[i]+0.001)
  }
  
  #------ correct names
  colnames(dl_max_min)[2:ncol(dl_max_min)] <- paste0(colnames(x)[2:ncol(x)], "_max_min")
  colnames(dl_log)[2:ncol(dl_log)] <- paste0(colnames(x)[2:ncol(x)], "_log")
  colnames(dl_rank)[2:ncol(dl_rank)] <- paste0(colnames(x)[2:ncol(x)], "_rank")

  
  #------ add total exposure for each normalization
  dl_max_min$expo_max_min <- dl_max_min$SR_tot_ias_max_min + 
    dl_max_min$med_range_max_min + dl_max_min$med_ib_max_min
  dl_log$expo_log <- dl_log$SR_tot_ias_log + 
    dl_log$med_range_log + dl_log$med_ib_log
  dl_rank$expo_rank <- dl_rank$SR_tot_ias_rank + 
    dl_rank$med_range_rank + dl_rank$med_ib_rank
  
  #------ test exposure as a product of components (for each normalization)
  dl_max_min$expo_prod_max_min <- dl_max_min$SR_tot_ias_max_min * 
    dl_max_min$med_range_max_min * dl_max_min$med_ib_max_min
  dl_log$expo_prod_log <- dl_log$SR_tot_ias_log * 
    dl_log$med_range_log * dl_log$med_ib_log
  dl_rank$expo_prod_rank <- dl_rank$SR_tot_ias_rank * 
    dl_rank$med_range_rank * dl_rank$med_ib_rank
  
  # bind all tables
  all_norm <- list(
    dl_max_min, dl_log, dl_rank) %>% 
    reduce(left_join, by="grid_id")
  
  return(all_norm)
  
})

lapply(dl_all_norm, summary)

hist(dl_all_norm$all_groups$SR_tot_ias_max_min)
hist(dl_all_norm$all_groups$SR_tot_ias_log)
hist(dl_all_norm$all_groups$SR_tot_ias_rank)

hist(dl_all_norm$all_groups$expo_max_min)
hist(dl_all_norm$all_groups$expo_log)
hist(dl_all_norm$all_groups$expo_rank)

hist(dl_all_norm$all_groups$expo_prod_max_min)
hist(dl_all_norm$all_groups$expo_prod_log)
hist(dl_all_norm$all_groups$expo_prod_rank)

# save exposure table
saveRDS(dl_all_norm, paste0("Output/Exposure/BIVA_10_expo_norm_", res, "_km"))





