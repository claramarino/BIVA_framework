# Compute sensitivity to IAS

rm(list=ls())

library(plyr)
library(tidyverse)
library(stringr)
library(readr)

# different possible resolutions (computed for 0.1 & 1 degree)
deg = "01" # chose 1 or 01

#### Load total SR ####

sr_fold <- "Output/Sensitivity/SR_Tot_native_all_BMR/"

list.files(sr_fold)

sr_b1 <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_BIRD_r", deg, "_chunk1"))
sr_b2 <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_BIRD_r", deg, "_chunk2"))
sr_b3 <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_BIRD_r", deg, "_chunk3"))
sr_b4 <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_BIRD_r", deg, "_chunk4"))
sr_m <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_MAM_r", deg))
sr_r <- readRDS(paste0(sr_fold, "RISK_32_SR_tot_REPT_r", deg))

sr_b <- bind_rows(sr_b1, sr_b2, sr_b3, sr_b4) %>%
  group_by(x, y) %>%
  summarise(SR_tot = sum (SR_tot))
rm(sr_b1, sr_b2, sr_b3, sr_b4)


sr_bmr <- bind_rows(sr_b, sr_m, sr_r) %>%
  group_by(x, y) %>%
  summarise(SR_tot_bmr = sum (SR_tot))

#### Load cells IAS-A ####

cell_fold <- "Output/Sensitivity/Cells_native_IAS_A/"
list.files(cell_fold)
cells_b <- readRDS(paste0(cell_fold, "RISK_31_cells_nat_ias_a_BIRD_r", deg))
cells_m <- readRDS(paste0(cell_fold, "RISK_31_cells_nat_ias_a_MAM_r", deg))
cells_r <- readRDS(paste0(cell_fold, "RISK_31_cells_nat_ias_a_REPT_r", deg))

df_iasa_b <- bind_rows(cells_b)
df_iasa_m <- bind_rows(cells_m)
df_iasa_r <- bind_rows(cells_r)

#### Get ias A/T info ####

# load threats
native_iasa_x_ias <- readRDS("characterize_ias_yan/Data/00_all_native_ias_a")
length(unique(native_iasa_x_ias$scientificName))
table(native_iasa_x_ias %>% distinct(scientificName, className) %>% pull(className))

# compute T, NT, EX-EW col
native_iasa <- native_iasa_x_ias %>%
  distinct(scientificName, className, redlistCategory) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Endangered", "Critically Endangered", "Vulnerable"),
    "T", "NT_DD")) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Extinct", "Extinct in the Wild"), 
    "EX_EW", threatened))


# add a column for species in ias_x_native with specified ias
native_spe_ias <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok") %>%
  distinct(scientificName) %>% pull(scientificName)

native_iasa <- native_iasa %>%
  mutate(ias_modeled = if_else(scientificName %in% native_spe_ias, "YES","NO"))


#### Compute metrics ####

# IAS-A all
sr_iasa_b <- df_iasa_b %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_b = n())
sr_iasa_m <- df_iasa_m %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_m = n())
sr_iasa_r <- df_iasa_r %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_r = n())

# IAS-A with IAS modeled 
sr_iasa_mod_b <- df_iasa_b %>%
  filter(binomial %in% native_spe_ias) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_mod_b = n())
sr_iasa_mod_m <- df_iasa_m %>%
  filter(binomial %in% native_spe_ias) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_mod_m = n())
sr_iasa_mod_r <- df_iasa_r %>%
  filter(binomial %in% native_spe_ias) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iasa_mod_r = n())

# IAS-T all
ias_t_sp <- native_iasa %>% filter(threatened == "T") %>% pull(scientificName)

sr_iast_b <- df_iasa_b %>%
  filter(binomial %in% ias_t_sp) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_b = n())
sr_iast_m <- df_iasa_m %>%
  filter(binomial %in% ias_t_sp) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_m = n())
sr_iast_r <- df_iasa_r %>%
  filter(binomial %in% ias_t_sp) %>% distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_r = n())

# IAS-T with IAS modeled 
sr_iast_mod_b <- df_iasa_b %>%
  filter(binomial %in% native_spe_ias & binomial %in% ias_t_sp) %>% 
  distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_mod_b = n())
sr_iast_mod_m <- df_iasa_m %>%
  filter(binomial %in% native_spe_ias & binomial %in% ias_t_sp) %>% 
  distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_mod_m = n())
sr_iast_mod_r <- df_iasa_r %>%
  filter(binomial %in% native_spe_ias & binomial %in% ias_t_sp) %>% 
  distinct() %>%
  group_by(x, y) %>% summarise(SR_iast_mod_r = n())

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

library(ggcorrplot)

mcor <- cor(sensit_tbl)

ggcorrplot(mcor[-c(1,2), -c(1,2)])

norm_sens <- sensit_tbl %>% 
  mutate(ias_a_norm = SR_iasa_m/SR_tot_m,
         ias_t_norm = SR_iast_m/SR_tot_m)
  

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


