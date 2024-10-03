# Calculate completeness components for sensitivity

rm(list=ls())

library(tidyverse)

# choose resolution 
res = "110" # 55 or 110 km

# open cell ids & sp presence per cell
sp_fold = "Output/Sensitivity/Cells_nat_all_bmr/"
cell_ids_m <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_MAM_all_nat_", res))
cell_ids_b1 <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_BIRD1_all_nat_", res))
cell_ids_b2 <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_BIRD2_all_nat_", res))
cell_ids_b3 <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_BIRD3_all_nat_", res))
cell_ids_b4 <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_BIRD4_all_nat_", res))
cell_ids_r <- readRDS(paste0(sp_fold, "RISK_33_cell_IDs_REPT_all_nat_", res))


cell_ids <- list(
  mam = cell_ids_m,
  bird = bind_rows(cell_ids_b1, cell_ids_b2, cell_ids_b3, cell_ids_b4) %>% distinct(),
  rept = cell_ids_r
)


# List all DD tetrapods
redl <- read.csv("characterize_ias_yan/Data/IUCN_download_23062022/simple_summary.csv")

table(redl$className, redl$redlistCategory)

# load associated natives
ias_x_native <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok")
# Load occurrence ias
in_fold = "Output/Exposure/Exposure_raw/"
all_sp_cells <- readRDS(paste0(in_fold, "RISK_21_grid_cells_310_IAS_", res, "km"))

table(ias_x_native$scope, ias_x_native$severity)
table(ias_x_native %>%
        # we keep only natives associated to the modeled ias 
        filter(ias_gbif_key %in% unique(all_sp_cells$new_key))%>%
        rename(new_key = ias_gbif_key) %>%
        filter(!redlistCategory %in% c("Extinct","Extinct in the Wild")) %>%
        distinct(scientificName, className) %>%
        pull(className))


ias_x_nat_dd_ss <- ias_x_native %>%
  # we keep only natives associated to the modeled ias 
  filter(ias_gbif_key %in% unique(all_sp_cells$new_key))%>%
  rename(new_key = ias_gbif_key) %>%
  filter(!redlistCategory %in% c("Extinct","Extinct in the Wild")) %>%
  mutate(sco = if_else(scope %in% c("", "Unknown"), 0, 1)) %>%
  mutate(sev = if_else(severity %in% c("", "Unknown"), 0, 1)) %>%
  group_by(scientificName, className) %>%
  summarise(sum_sco = sum(sco),
            sum_sev = sum(sev)) %>% 
  filter(sum_sco<1 & sum_sev <1)
# all species badly informed for the IAS threat (severity /scope / ias name)
table(ias_x_nat_dd_ss$className)
# 18 birds, 166 mammals, 203 reptiles

dd_ias_a <- ias_x_nat_dd_ss %>% pull(scientificName)

# chargement de tous les IUCN DD = 
# toutes les sp qui sont mal renseignées en général dans l'iucn
fold_iucn <- "Z:/THESE/5_Data/IUCN_summary/MAJ_2023_06/" 
bp <- read.csv(paste0(fold_iucn, "birds_passerif/simple_summary.csv"))
bnp <- read.csv(paste0(fold_iucn, "birds_nonpasserif/simple_summary.csv"))
amr <- read.csv(paste0(fold_iucn, "amph_mam_rept/simple_summary.csv"))

all <- bind_rows(bp, bnp, amr) %>%
  filter(className != "AMPHIBIA")
table(all$redlistCategory, all$className)
table(all$className)

dd_bmr <- all %>% filter(redlistCategory=="Data Deficient") %>% pull(scientificName)

# 46 birds, 839 mammals, 1487 reptiles
table(all %>% distinct(scientificName, className) %>% pull(className))

# calculate SR for each condition

sr_iasa <- readRDS(paste0("Output/Sensitivity/SR_per_cell/RISK_33_SR_tot_ias_a_t_nt_", res))

sr_dd_iasa <- lapply(cell_ids, function(x){
  x %>% 
    filter(binomial %in% dd_ias_a) %>%
    group_by(cell_id) %>%
    summarise(SR_dd_iasa = n())
})

sr_dd_all <- lapply(cell_ids, function(x){
  x %>% 
    filter(binomial %in% dd_bmr) %>%
    group_by(cell_id) %>%
    summarise(SR_dd_all = n())
})

sr_all <- mapply(left_join, 
                 mapply(full_join, sr_iasa, sr_dd_iasa, SIMPLIFY = FALSE), 
                 sr_dd_all, 
                 SIMPLIFY = FALSE)
colnames(sr_all$mam)

sr_all <- lapply(sr_all, function(x){
  x %>% replace_na(list(SR_dd_all = 0, SR_dd_iasa = 0)) %>%
    filter(!is.na(cell_id))
})
lapply(sr_all, cor)


saveRDS(sr_all, paste0("Output/Sensitivity/Completeness/RISK_34_SR_DD_sensit_", res))


###### Calculate completeness 

sr_all <- readRDS(paste0("Output/Sensitivity/Completeness/RISK_34_SR_DD_sensit_", res))

sr_all$all <- bind_rows(sr_all) %>% group_by(cell_id) %>%
  summarise_at(c("SR_tot", "SR_ias_a", "SR_dd_iasa", "SR_dd_all"), sum, na.rm = TRUE)

compl <- lapply(sr_all, function (x){x %>%
    mutate(se = 1-(SR_dd_all/SR_tot), 
           ed = 1-(SR_dd_iasa/SR_ias_a)) %>%
    mutate(comp_prod = se * ed)})


hist(compl$all$se)
hist(compl$all$ed)
hist(compl$all$comp_prod)

hist(compl$rept$se)
hist(compl$rept$ed)
hist(compl$rept$comp_prod)

hist(compl$mam$se)
hist(compl$mam$ed)
hist(compl$mam$comp_prod)

hist(compl$bird$se)
hist(compl$bird$ed)
hist(compl$bird$comp_prod)

saveRDS(compl, paste0("Output/Sensitivity/Completeness/RISK_34_completeness_ed_se_", res))
