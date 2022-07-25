# test if exposure to direct, indirect or habitat are congruent

rm(list=ls())

library(tidyverse)
library(sp)
library(sf)
library(raster)


# Load exposure data
out_fold = "Output/Raster_metrics/Exposure_raw/"
all_sp_cells <- readRDS(paste0(out_fold, "RISK_21_raster_cells_308_IAS_r1"))
# combine df for each sp into one df
df_all <- bind_rows(all_sp_cells) 

# for each ias, calculate max range
tot_range <- df_all %>%
  group_by(new_key, new_species) %>%
  summarise(n_cells_range = n())

df_all_range <- left_join(df_all %>% dplyr::select(-new_key), 
                          tot_range,
                          by = "new_species")

# load ias info 
ias_data <- readRDS("Output/IAS_characteristics/RISK_21_IAS_meca_asso_nat_df")

quantile(ias_data$IAS_A_tot)
quantile(tot_range$n_cells_range, c(0.10,0.9))
hist((tot_range$n_cells_range), breaks =50)
abline(v=592)


ias_data_summary <- ias_data %>% 
  mutate()

df_all_data <- left_join(
  df_all_range %>% rename(ias_gbif_key = new_key),
  ias_data_summary, 
  by = "ias_gbif_key")




# calculate cell related metrics
ias_all_agg <- df_all_range %>%
  group_by(x, y, cells) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_pts), # nb occ per cell (all ias sp)
    nb_occ_med = median(n_pts), # median nb occ per cell
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range)
  ) 

coast <- rnaturalearth::ne_coastline(returnclass = "sf")

# plot species richness
SR_ias<- ggplot(data = ias_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  geom_sf(data = coast, alpha = 0.1) +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$new_key))),
       x = "Longitude", y = "Latitude")
SR_ias


#### calculate metrics per meca ####

meca_agg <- df_all_data %>%
  group_by(x, y, cells, major_meca) %>%
  summarise(
    SR_tot = n(), # alien species richness per cell
    nb_occ_tot = sum(n_pts), # nb occ per cell (all ias sp)
    nb_occ_med = median(n_pts), # median nb occ per cell
    range_tot = sum(n_cells_range),
    range_med = median(n_cells_range)) 

meca_agg_wide <- meca_agg %>%
  dplyr::select(x:SR_tot) %>%
  filter(!is.na(major_meca)) %>%
  pivot_wider(names_from = major_meca, values_from = SR_tot, values_fill = 0) %>%
  mutate(meca_cell = case_when(
    Ecosystem > 0 & Indirect_sp_effect == 0 & Direct_sp_effect == 0 ~ "Only_eco",
    Ecosystem == 0 & Indirect_sp_effect > 0 & Direct_sp_effect == 0 ~ "Only_indir",
    Ecosystem == 0 & Indirect_sp_effect == 0 & Direct_sp_effect > 0 ~ "Only_dir",
    Ecosystem > 0 & Indirect_sp_effect > 0 & Direct_sp_effect == 0 ~ "Eco_indir",
    Ecosystem > 0 & Indirect_sp_effect == 0 & Direct_sp_effect > 0 ~ "Eco_dir",
    Ecosystem == 0 & Indirect_sp_effect > 0 & Direct_sp_effect > 0 ~ "Dir_indir",
    Ecosystem > 0 & Indirect_sp_effect > 0 & Direct_sp_effect > 0 ~ "all"
  ))

meca_agg_wide
meca<- ggplot(data = meca_agg_wide) +
  geom_raster(aes(x = x, y = y, fill = meca_cell)) +
  geom_sf(data = coast, alpha = 0.1) +
  theme_classic() +
  labs(title = "Dominant meca within each cell",
       x = "Longitude", y = "Latitude")
meca

cor_meca_ei <- ggplot(data = meca_agg_wide, 
                   aes(x = Ecosystem, y = Indirect_sp_effect)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope=1, intercept = 0)
cor_meca_ei             

cor_meca_ed <- ggplot(data = meca_agg_wide, 
                   aes(x = Ecosystem, y = Direct_sp_effect)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope=1, intercept = 0)
cor_meca_ed   

cor_meca_id <- ggplot(data = meca_agg_wide, 
                      aes(x = Direct_sp_effect, y = Indirect_sp_effect)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope=1, intercept = 0)
cor_meca_id   


table(ias_data %>% filter(ias_gbif_key %in% unique(df_all$new_key)) %>% 
        pull(major_meca))



##### Calculate metrics using IAS-T, IAs-EX, IAS-NT #######
