# compute exposure metrics

# one metric with species richness
# one metric accounting for ias range
# one metric accounting for assoc natives

# normalize
# save table for analysis with sensitivity

rm(list=ls())

library(tidyverse)
library(sp)
library(sf)
library(raster)


# select resolution
deg_in = "0.1" # "1" or "01"
deg_out = "01"

# Load data
in_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(in_fold, "RISK_21_raster_cells_308_IAS_r", deg_in))

# combine df for each sp into one df
df_all <- bind_rows(all_sp_cells) 

str(df_all)

# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())

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
  group_by(x, y, cells) %>%
  summarise(
    # metric 1 = alien species richness per cell
    SR_tot_ias = n(), 
    # metric 2 = account for species range
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range),
    range_mean = mean(n_cells_range),
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


mcor <- cor(na.omit(exposure %>% dplyr::select(-c(cells, mean_sco_sev))))
ggcorrplot(mcor)


exposure_tbl <- exposure %>%
  dplyr::select(x,y, SR_tot_ias, range_med, med_iasa_tot)


# normalize exposure metrics

# max min linear
maxcol <- apply(exposure_tbl, 2, max)
mincol <- apply(exposure_tbl, 2, min)

expo_tbl_norm_lin <- exposure_tbl
for (i in 3:length(maxcol)){
  expo_tbl_norm_lin[,i] <- (expo_tbl_norm_lin[,i]-mincol[i])/(maxcol[i]-mincol[i])
}
summary(expo_tbl_norm_lin)


#log transform and max min linear
maxlog <- apply(exposure_tbl, 2, function(x){max(log(x+1))})
minlog <- apply(exposure_tbl, 2, function(x){min(log(x+1))})
expo_tbl_norm_log <- exposure_tbl
for (i in 3:length(maxlog)){
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

for (i in 3:length(maxcol)){
  expo_tbl_rank[,i] <- (expo_tbl_rank[,i]-minrank[i])/(maxrank[i]-minrank[i])
}

summary(expo_tbl_rank)



# save exposure tables
saveRDS(expo_tbl_norm_lin, paste0("Output/Exposure/RISK_24_expo_norm_lin_r", deg_out))
saveRDS(expo_tbl_norm_log, paste0("Output/Exposure/RISK_24_expo_norm_log_r", deg_out))
saveRDS(expo_tbl_rank, paste0("Output/Exposure/RISK_24_expo_norm_rank_r", deg_out))



