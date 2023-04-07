# compute exposure metrics
# following meeting with guillaume

# normalize
# save table for analysis with sensitivity

rm(list=ls())

library(tidyverse)
library(sp)
library(sf)


###### Data loading #####

# select resolution
nkm = "55km" # 55km or 110km

# Load occurrence ias
in_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(in_fold, "RISK_21_grid_cells_310_IAS_", nkm))
# cb de cellules avec au moins 1 IAS ?
length(unique(all_sp_cells$grid_id))

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
  # filter(!redlistCategory %in% c("Extinct","Extinct in the Wild")) %>%
  dplyr::select(new_key, scientificName, className, redlistCategory) %>%
  distinct()

sum(unique(all_sp_cells$new_key) %in% unique(ias_x_native_ok$new_key))

###### Compute cell related metrics #####

# combine for calculating cell related metrics
df_all <- left_join(all_sp_cells,
                    ias_x_native_ok,
                    relationship = "many-to-many")

str(df_all)


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

# impact breadth metric 1
# function for counting associated victim SR per cell,
# with possible redundancy if victim is associated to several ias present
# union + intersection of impact breadth
calc_ib1_red <- function(x){
  out <- x %>%
    dplyr::select(grid_id, scientificName) %>%
    # distinct() %>% # no distinct bc we want possible duplicates
    group_by(grid_id) %>%
    summarise(ib1_red = n())
  return(out)
}

# impact breadth metric 2
# function for counting associated victim SR per cell,
# without redundancy = union of impact breadth
calc_ib2_no_red <- function(x){
  out <- x %>%
    dplyr::select(grid_id, scientificName) %>%
    distinct() %>%
    group_by(grid_id) %>%
    summarise(ib2_no_red = n())
  return(out)
}

# impact union metric 
# function for counting the victim sp that are associated to several IAS
# remove lonely species bc they are already in IB
calc_iu <- function(x){
  temp <- x %>%
    group_by(grid_id, scientificName) %>%
    summarise(dupl_vict = n())
  out <- temp %>%
    filter(dupl_vict>1) %>%
    group_by(grid_id) %>%
    summarise(iu = sum(dupl_vict))
  return(out)
}


# apply functions to the list with all taxa and 1 df for each victim group

dl_metrics <- lapply(dl, function (df){
  met_to_join <- list(
    calc_sr_ias(df),
    calc_ib1_red(df),
    calc_ib2_no_red(df),
    calc_iu(df))
  out <- met_to_join %>% reduce(left_join, by="grid_id") %>%
    mutate(iu = replace_na(iu, 0))
  return(out)
})



library(ggcorrplot)
mcor <- cor(dl_metrics$all_groups %>% 
              dplyr::select(-grid_id)%>%
              dplyr::mutate(ibu = iu+ib2_no_red))
ggcorrplot(mcor)

hist(dl_metrics$all_groups %>%
       dplyr::mutate(ibu = iu+ib2_no_red) %>%
       pull(ibu))

hist(dl_metrics$all_groups %>%
       pull(iu))

hist(dl_metrics$all_groups %>%
       pull(ib1_red))

hist(dl_metrics$all_groups %>%
       pull(ib2_no_red))

hist(dl_metrics$all_groups %>%
       pull(SR_tot_ias))


ggplot(dl_metrics$rept, aes(x=SR_tot_ias, y = ib2_no_red))+
  geom_point()

ggplot(dl_metrics$rept, aes(x=SR_tot_ias, y = ib1_red))+
  geom_point()+
  geom_smooth()

ggplot(dl_metrics$rept, aes(x=SR_tot_ias, y = iu))+
  geom_point()+
  geom_smooth()

nb_nat_ias <- ias_x_native_ok %>% count(new_key)
table(ias_x_native_ok$redlistCategory)

####################### BROUILLON ########################

# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n(),
            area_range = sum(area)) %>%
  mutate(ln_n_cells = log(n_cells_range),
         ln_area_range = log(area_range))

df_all_range <- left_join(df_all %>% dplyr::select(-new_key), 
                          tot_range,
                          by = "new_species")


# load IAS info on severity and IAS-A
ias_info <- readRDS("Output/IAS_characteristics/RISK_22_IAS_meca_sev_asso_nat_df")

df_all_info <- left_join(df_all_range, 
                         ias_info %>% dplyr::select(-ias_species) %>%
                           rename(new_key = ias_gbif_key))


# calculate cell related metrics
exposure <- df_all_info %>%
  group_by(x, y, cells, area) %>%
  summarise(
    # metric 1 = alien species richness per cell
    SR_tot_ias = n(), 
    # metric 2 = account for species range
    range_tot = sum(ln_area_range),
    range_med = median(ln_area_range),
    range_mean = mean(ln_area_range),
    # metric 3 = account for IAS-A, T, EX
    sum_iasa_tot = sum(IAS_A_tot),
    sum_iast = sum(`T`),
    sum_iasnt = sum(NT_DD),
    sum_iasex = sum(EX_EW),
    sum_sco_sev = sum(sco_sev_ias, na.rm = T),
    med_iasa_tot = median(IAS_A_tot),
    med_iast = median(`T`),
    med_iasnt = median(NT_DD),
    med_iasex = median(EX_EW),
    mean_iasa_tot = mean(IAS_A_tot),
    mean_iast = mean(`T`),
    mean_iasnt = mean(NT_DD),
    mean_iasex = mean(EX_EW),
    mean_sco_sev = mean(sco_sev_ias, na.rm = T)
  ) 

library(ggcorrplot)
mcor <- cor(na.omit(exposure %>% dplyr::select(-c(cells, mean_sco_sev))))
ggcorrplot(mcor)


exposure_tbl <- exposure %>%
  dplyr::select(x, y, area, SR_tot_ias, range_med, med_iasa_tot) %>%
  mutate(SR_tot_ias_area = (SR_tot_ias/area)*1000)
# multiply by 1000 to smooth the repartition (because dividing by area => small values)

cor.test(exposure_tbl$SR_tot_ias, exposure_tbl$SR_tot_ias_area)
hist(exposure_tbl$SR_tot_ias_area)

# normalize exposure metrics

# max min linear
maxcol <- apply(exposure_tbl, 2, max)
mincol <- apply(exposure_tbl, 2, min)

expo_tbl_norm_lin <- exposure_tbl
for (i in 5:length(maxcol)){
  expo_tbl_norm_lin[,i] <- (expo_tbl_norm_lin[,i]-mincol[i])/(maxcol[i]-mincol[i])
}
summary(expo_tbl_norm_lin)


#log transform and max min linear
maxlog <- apply(exposure_tbl, 2, function(x){max(log(x+1))})
minlog <- apply(exposure_tbl, 2, function(x){min(log(x+1))})
expo_tbl_norm_log <- exposure_tbl
for (i in 5:length(maxlog)){
  expo_tbl_norm_log[,i] <- (log(expo_tbl_norm_log[,i]+1)-minlog[i])/
    (maxlog[i]-minlog[i])
}
summary(expo_tbl_norm_log)


# cumulative distribution
expo_tbl_rank <- exposure_tbl %>%
  dplyr::ungroup() %>%
  dplyr::mutate_at(vars(-c(x,y)), dense_rank)
maxrank <- apply(expo_tbl_rank, 2, max)
minrank <- apply(expo_tbl_rank, 2, min)

for (i in 5:length(maxcol)){
  expo_tbl_rank[,i] <- (expo_tbl_rank[,i]-minrank[i])/(maxrank[i]-minrank[i])
}

summary(expo_tbl_rank)

hist(expo_tbl_norm_lin$SR_tot_ias_area)
hist(expo_tbl_norm_log$SR_tot_ias_area)
hist(expo_tbl_rank$SR_tot_ias_area)


# save exposure tables
saveRDS(expo_tbl_norm_lin, paste0("Output/Exposure/RISK_24_expo_norm_lin_r", deg_out))
saveRDS(expo_tbl_norm_log, paste0("Output/Exposure/RISK_24_expo_norm_log_r", deg_out))
saveRDS(expo_tbl_rank, paste0("Output/Exposure/RISK_24_expo_norm_rank_r", deg_out))



