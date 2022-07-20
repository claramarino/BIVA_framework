# count natives associated to IAS
# find associated mecha to each IAS

rm(list=ls())

library(tidyverse)

# load ias list
sp_info <- readRDS("Output/Native_exotic_range/RISK_14_ias_list_with_occ_ALL_DB_nat_exo")

# load ias_x_native table
ias_x_native <- readRDS("Data/RISK_03_ias_x_native_bmr_taxo_ok")


length(unique(ias_x_native$ias))
length(unique(ias_x_native$ias_gbif_key))


##### Number of IAS-A / IAS-T associated to each IAS #####

ias_nat_count <- ias_x_native %>%
  select(className, scientificName, redlistCategory,
         ias_gbif_key, ias_species) %>% 
  distinct() %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Endangered", "Critically Endangered", "Vulnerable"),
    "T", "NT_DD")) %>%
  mutate(threatened = if_else(
    redlistCategory %in% c("Extinct", "Extinct in the Wild"), 
    "EX_EW", threatened)) %>%
  count(ias_gbif_key, ias_species, threatened) %>%
  pivot_wider(names_from = threatened, values_from = n, values_fill = 0) %>%
  mutate(IAS_A_tot = `T` + NT_DD + EX_EW)


##### IAS mechanisms #####

# remove column with wrong ias_names (before taxo check)
# keep only native info + ias key + meca
# deal with meca column
ias_x_native_meca <- ias_x_native %>%
  select(className, scientificName, redlistCategory,
         ias_gbif_key, ias_species, stressName) %>% 
  distinct() %>%
  separate(stressName, c("nMeca1","nMeca2","nMeca3","nMeca4", "nMeca5"), '([|])', 
           extra = 'merge') %>%
  pivot_longer(cols = starts_with("nMeca"),
               names_to = "Meca_num",
               values_to = "Meca_name",
               values_drop_na = TRUE) %>%
  # remove empty mechanism
  filter(!Meca_name %in% c("", "Other"))%>%
  mutate(Cate_mecha = case_when(
    Meca_name=="Ecosystem conversion" | Meca_name=="Ecosystem degradation" |
      Meca_name=="Indirect ecosystem effects" ~ "Ecosystem",
    Meca_name=="Species mortality" |  Meca_name=="Reduced reproductive success"|
      Meca_name=="Species disturbance" ~  "Direct_sp_effect",
    Meca_name=="Hybridisation"| Meca_name=="Competition" |Meca_name=="Inbreeding" |
      Meca_name=="Skewed sex ratios"~ "Indirect_sp_effect")) %>%
  select(-c(Meca_num, Meca_name)) %>% distinct()

ias_mecha_count <- ias_x_native_meca %>%
  count(ias_gbif_key, Cate_mecha) %>%
  pivot_wider(names_from = Cate_mecha, values_from = n, values_fill = 0) 
ias_mecha_count_maj <- ias_mecha_count %>%
  mutate(major_meca = colnames(ias_mecha_count %>% 
                                 select(Ecosystem:Indirect_sp_effect)
                               )[max.col(ias_mecha_count %>%
                                           select(Ecosystem:Indirect_sp_effect),
                                         ties.method = "first")])




# final IAS info 

ias_data_nat_mec <- left_join(ias_nat_count, ias_mecha_count_maj, 
                              by ="ias_gbif_key") %>%
  # arrange columns
  select(ias_gbif_key, ias_species, IAS_A_tot, NT_DD:EX_EW, major_meca, Ecosystem:Indirect_sp_effect)

saveRDS(ias_data_nat_mec, "Output/IAS_characteristics/RISK_21_IAS_meca_asso_nat_df")
