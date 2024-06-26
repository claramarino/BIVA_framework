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

# add number of IAS-AS

ias_as <- ias_x_native %>% 
  distinct(scientificName, ias_gbif_key, ias_species, scope, severity) %>%
  mutate(IAS_AS = if_else(
    scope %in% c("Majority (50-90%)","Whole (>90%)") & 
      severity %in% c("Slow, Significant Declines", "Rapid Declines", "Very Rapid Declines"), 
    "IAS_AS", "no")) %>%
  count(ias_gbif_key, ias_species, IAS_AS) %>%
  pivot_wider(names_from = IAS_AS, values_from = n, values_fill = 0) %>%
  select(ias_gbif_key, ias_species, IAS_AS)

# 125 IAS with 1 or + native severely associated

ias_nat_count <- left_join(ias_nat_count, ias_as)

##### CHECK per species effect and EICAT impact ####

# EICAT
eicat <- read.csv2("Z:/THESE/6_Projects/mechanisms_birds/Data/EICAT_birds_2016.csv") %>%
  mutate(binomial = paste(Gender, Species, sep=" ")) %>%
  select(binomial, Impact_cate, Impact_mecha, Assessment)

eicat_nat <- left_join(eicat, ias_nat_count %>% rename(binomial = ias_species),
                       by = "binomial") %>%
  filter(!is.na(ias_gbif_key))

ggplot(data = eicat_nat, aes(x = Impact_cate, y = IAS_A_tot))+
  geom_point(position = "jitter")
ggplot(data = eicat_nat) +
  # geom_point(aes(x = Impact_cate, y = IAS_A_tot),
  #            color = "black", position = "jitter") +
  geom_point(aes(x = Impact_cate, y = NT_DD),
             color = "turquoise", position = "jitter", show.legend = T) +
  geom_point(aes(x = Impact_cate, y = `T`),
             color = "orange", position = "jitter", show.legend = T) +
  geom_point(aes(x = Impact_cate, y = EX_EW),
             color = "firebrick", position = "jitter", show.legend = T)


eicat_nat_l <- eicat_nat %>% pivot_longer(
  cols = NT_DD:IAS_AS,
  names_to = "group_nat",
  values_to = "nb_sp") %>%
  mutate(Impact_cate2 = if_else(Impact_cate %in% c('MO','MR','MV'), "harmful", "_non_harmful")) %>%
  mutate(Impact_cate2 = if_else(Impact_cate =="DD", "_DD", Impact_cate2))

ggplot(data = eicat_nat_l %>% 
         filter(group_nat!="IAS_A_tot") %>%
         filter(nb_sp>0))+
  geom_boxplot(aes(x = Impact_cate2, y = nb_sp, color = group_nat))
  
ggplot(data = eicat_nat_l %>% 
         filter(group_nat!="IAS_A_tot") %>%
         filter(nb_sp>0))+
  geom_point(aes(x = Impact_cate2, y = nb_sp, color = group_nat), alpha=0.5, position ="jitter")


ggplot(data = eicat_nat_l %>% 
         filter(group_nat=="IAS_A_tot"))+
  geom_boxplot(aes(x = Impact_cate, y = nb_sp))+
  geom_point(aes(x = Impact_cate, y = nb_sp), alpha=0.5, position ="jitter")


ggplot(data = eicat_nat_l %>% 
         filter(group_nat=="IAS_AS"))+
  geom_boxplot(aes(x = Impact_cate, y = nb_sp))+
  geom_point(aes(x = Impact_cate, y = nb_sp), alpha=0.5, position ="jitter")

ggplot(data = eicat_nat_l %>% 
         filter(group_nat=="IAS_AS") %>%
         filter(nb_sp>0))+
  geom_boxplot(aes(x = Impact_cate, y = nb_sp))+
  geom_point(aes(x = Impact_cate, y = nb_sp), alpha=0.5, position ="jitter")

ggplot(data = eicat_nat_l %>% 
         filter(group_nat=="IAS_AS") %>%
         filter(nb_sp>0))+
  geom_boxplot(aes(x = Impact_cate2, y = nb_sp))+
  geom_point(aes(x = Impact_cate2, y = nb_sp), alpha=0.5, position ="jitter")

ggplot(data = eicat_nat_l %>% 
         filter(group_nat=="IAS_AS") %>%
         filter(nb_sp>0))+
  geom_boxplot(aes(x = Impact_cate, y = nb_sp))+
  geom_point(aes(x = Impact_cate, y = nb_sp), alpha=0.5, position ="jitter")

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



###### Associations meca x extinction ######

asso_meca_threat <- ias_x_native_meca %>% 
  select(className, scientificName, redlistCategory, Cate_mecha) %>%
  mutate(count = 1) %>%
  pivot_wider(names_from = Cate_mecha,
              values_from = count,
              values_fn = sum, 
              values_fill = 0) %>%
  mutate(combi = case_when(
    Ecosystem > 0 & Indirect_sp_effect == 0 & Direct_sp_effect == 0 ~ "Only_eco",
    Ecosystem == 0 & Indirect_sp_effect > 0 & Direct_sp_effect == 0 ~ "Only_indir",
    Ecosystem == 0 & Indirect_sp_effect == 0 & Direct_sp_effect > 0 ~ "Only_dir",
    Ecosystem > 0 & Indirect_sp_effect > 0 & Direct_sp_effect == 0 ~ "Eco_indir",
    Ecosystem > 0 & Indirect_sp_effect == 0 & Direct_sp_effect > 0 ~ "Eco_dir",
    Ecosystem == 0 & Indirect_sp_effect > 0 & Direct_sp_effect > 0 ~ "Dir_indir",
    Ecosystem > 0 & Indirect_sp_effect > 0 & Direct_sp_effect > 0 ~ "all"
  )) %>%
  filter(redlistCategory!="Data Deficient") %>%
  mutate(ex_risk = case_when(
    redlistCategory == "Least Concern" ~ 1, 
    redlistCategory == "Near Threatened" ~ 2,
    redlistCategory == "Vulnerable" ~ 3,
    redlistCategory == "Endangered" ~ 4,
    redlistCategory == "Critically Endangered" ~ 5,
    redlistCategory == "Extinct in the Wild" ~ 6,
    redlistCategory == "Extinct" ~ 7))

table(asso_meca_threat$combi, asso_meca_threat$ex_risk)

ggplot(asso_meca_threat) +
  geom_boxplot(aes(x = combi, y = ex_risk)) +
  geom_point(aes(x = combi, y = ex_risk), position = "jitter")


asso_meca_threat %>% group_by(combi) %>% 
  summarise(mean_risk = mean(ex_risk))

m_class <- asso_meca_threat %>% group_by(combi, className) %>% 
  summarise(mean_risk = mean(ex_risk))

table(ias_x_native_meca$redlistCategory)


### add scope severity

scope_sev <- ias_x_native %>% 
  distinct(scientificName, ias_gbif_key, scope, severity) %>%
  filter(!(scope == "" & severity == "")) %>%
  filter(!(scope == "Unknown" & severity == "Unknown")) %>%
  mutate(scope_num = case_when(
    scope == "Unknown" | scope == "" ~ 0,
    scope == "Minority (<50%)" ~ 1,
    scope == "Majority (50-90%)" ~ 2,
    scope == "Whole (>90%)" ~ 3)) %>%
  mutate(sev_num = case_when(
    severity == "Unknown" | severity == "" ~ 0,
    severity == "No decline" ~ 1,
    severity == "Negligible declines" ~ 2,
    severity == "Causing/Could cause fluctuations" ~ 3,
    severity == "Slow, Significant Declines" ~ 4,
    severity == "Rapid Declines" ~ 5,
    severity == "Very Rapid Declines" ~ 6)) %>%
  mutate(scope_norm = scope_num/3,
         sev_norm = sev_num/6) %>%
  mutate(sco_sev = scope_norm + sev_norm)

head(scope_sev)

# compute one value of severity per IAS
ias_sco_sev <- scope_sev %>% group_by(ias_gbif_key) %>%
  summarize(sco_sev_ias = max(sco_sev),
            med_sco_sev = median(sco_sev),
            sco = max(scope_norm),
            sev = max(sev_norm))

###### Link with EICAT

eicat_nat <-  eicat_nat %>%
  mutate(Impact_cate2 = if_else(Impact_cate %in% c('MO','MR','MV'), "harmful", "_non_harmful")) %>%
  mutate(Impact_cate2 = if_else(Impact_cate =="DD", "_DD", Impact_cate2))

# check if ias scosev is consistent with EICAT score
ias_sev_eicat <- inner_join(eicat_nat, ias_sco_sev)
ggplot(data = ias_sev_eicat)+
  geom_boxplot(aes(x = Impact_cate, y = sco_sev_ias))+
  geom_point(aes(x = Impact_cate, y = sco_sev_ias), position = "jitter")

ggplot(data = ias_sev_eicat %>% filter(Assessment!="Low"))+
  geom_boxplot(aes(x = Impact_cate2, y = sco_sev_ias))+
  geom_point(aes(x = Impact_cate2, y = sco_sev_ias), position = "jitter")

ggplot(data = ias_sev_eicat)+
  geom_boxplot(aes(x = Impact_cate2, y = med_sco_sev))+
  geom_point(aes(x = Impact_cate2, y = med_sco_sev), position = "jitter")

ggplot(data = ias_sev_eicat)+
  geom_boxplot(aes(x = Impact_cate2, y = sco))+
  geom_point(aes(x = Impact_cate2, y = sco), position = "jitter")
ggplot(data = ias_sev_eicat)+
  geom_boxplot(aes(x = Impact_cate2, y = sev))+
  geom_point(aes(x = Impact_cate2, y = sev), position = "jitter")









# join with IAS meca
ias_sev_mecha <- inner_join(ias_mecha_count_maj, ias_sco_sev)
# link meca /severity
ggplot(ias_sev_mecha) +
  geom_boxplot(aes(major_meca, sco_sev_ias)) +
  geom_point(aes(major_meca, sco_sev_ias), position = "jitter")

one.way <- aov(sco_sev_ias ~ major_meca, data = ias_sev_mecha)

summary(one.way)
plot(one.way)

# compute one value of severity per native species
nat_sco_sev <- scope_sev %>% group_by(scientificName) %>%
  summarize(sco_sev_nat = mean(sco_sev),
            max_sco_sev = max(sco_sev))
# join with native meca
nat_sev_meca <- inner_join(asso_meca_threat, nat_sco_sev)
ggplot(nat_sev_meca) +
  geom_boxplot(aes(combi, max_sco_sev)) +
  geom_point(aes(combi, max_sco_sev), position = "jitter")

unique(scope_sev$severity)
table(scope_sev$severity)

table(asso_meca_threat$className)

# save severity scope info
ias_data_nat_mec_sev <- left_join(ias_data_nat_mec, ias_sco_sev,
                              by ="ias_gbif_key")
saveRDS(ias_data_nat_mec_sev, 
        "Output/IAS_characteristics/RISK_22_IAS_meca_sev_asso_nat_df")

summary(ias_data_nat_mec_sev)
