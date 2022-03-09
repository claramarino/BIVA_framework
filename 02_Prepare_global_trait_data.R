rm(list=ls())


# Load packages
# Automatically install required packages, which are not yet in library
packages <- c("tidyverse","BAT", "Ternary")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)


# load threat assemblages and IUCN categories
setwd("Z:/THESE/6_Projects/insular_tetrapods_FD")
df_all <- readRDS("Output/Data_clean/01_1_df_all_tetrapods_MAJ")
threat_ias <- readRDS("Output/Data_clean/01_1_threat_ias_all_tetrapods_MAJ") %>%
  rename(binomial= scientificName)
all_threats <- readRDS("Output/Data_clean/01_1_other_threat_all_tetra") %>%
  rename(binomial= scientificName)

# load native species traits for mammals & squamates (insular_tetrapods)
df_mam_5t_imputed <- readRDS("Output/Data_clean/02_3_1000_imp_mam_case2_15_10_MAJ")
df_rept_5t_imputed <- readRDS("Output/Data_clean/02_3_1000_imp_rept_case2_15_10_MAJ")

# load native traits for birds (mecha_birds)
setwd("Z:/THESE/6_Projects/mechanisms_birds")
df_birds <- readRDS("Output/Data_clean/08_df_all_imp_final")

df_birds_ft <- df_birds %>%
  mutate_if(is.character, as.factor) %>%
  column_to_rownames("binomial") %>%
  # convert Mass, Beak & Clutch with log 
  mutate(ln.Mass = log(Mass),
         ln.Clutch = log(Clutch),
         ln.Beak.Depth = log(Beak.Depth),
         ln.Beak.Length_Nares = log(Beak.Length_Nares)) %>%
  # remove converted var + trophic niche
  select(-c(Mass, Clutch, Beak.Depth, Beak.Length_Nares, Trophic.Niche))

# -------------------------------- MAMMALS 

df_out <- df_mam_5t_imputed[[1]]
for(i in 2:ncol(df_out)){
  df_out[,i] <- 
    apply(do.call("rbind",lapply(df_mam_5t_imputed,"[[", i)), 2, median)
} 

df_out$Diet_plant <- 0
df_out$Diet_vert <- 0
df_out$Diet_inv <- 0
for (i in 1:nrow(df_out)){
  df_out$Diet_plant[i] <- round(XYToTernary(
    df_out$Diet.x[i], df_out$Diet.y[i],1)[1]*100,0)
  df_out$Diet_vert[i] <- round(XYToTernary(
    df_out$Diet.x[i], df_out$Diet.y[i],1)[2]*100,0)
  df_out$Diet_inv[i] <- round(XYToTernary(
    df_out$Diet.x[i], df_out$Diet.y[i],1)[3]*100,0)
}

data_ft_mam <- df_out %>% 
  # remove species with missing values
  filter(!(is.na(Activity.Nocturnal) | 
             is.na(ln.Mass.g) | is.na(Diet.x) | `1:m_nbNA`==1))%>%
  # select only insular endemic species
  # filter(insular_endemic==1) %>% # keep all native mammals
  # Habitat trait
  mutate(hab_sum=as.character(hab_sum)) %>%
  mutate(hab_sum=if_else(hab_sum>4,"5",hab_sum)) %>%
  mutate(Habitat=factor(hab_sum, ordered = T,
                        levels = c("1","2","3","4","5"))) %>%
  # Body mass
  mutate(Body.mass=factor(ntile(ln.Mass.g,5), ordered = T,
                          levels = c("1","2","3","4","5"))) %>%
  # Diet 
  mutate(Diet.anim = Diet_inv + Diet_vert) %>%
  mutate(Main.diet = case_when(
    Diet_plant > 50 ~ "Herbi",
    Diet.anim > 50 ~ "Multi_anim",
    Diet.anim > 45 & Diet_plant > 45 ~ "Omni")) %>% # not ==50 because of imputation
  mutate(Main.diet=if_else(Diet_vert > 50, "Vert", Main.diet)) %>%
  mutate(Main.diet=if_else(Diet_inv > 50, "Inv", Main.diet)) %>%
  # Foraging niche
  mutate(For.niche = case_when(
    `1:m_nbAr`>0 ~ "Arboreal",
    `1:m_nbG`>0 ~ "Ground",
    `1:m_nbM`>0 ~ "Marine",
    `1:m_nbA`>0 ~ "Aerial",
    `1:m_nbS`>0 ~ "Scansorial")) %>%
  # foraging period (activity)
  mutate(Activity = case_when(
    Activity.Nocturnal==1 & Activity.Crepuscular==1 & Activity.Diurnal==1~ "all",
    Activity.Nocturnal==0 & Activity.Crepuscular==1 & Activity.Diurnal==1~ "CD",
    Activity.Nocturnal==1 & Activity.Crepuscular==1 & Activity.Diurnal==0~ "CN",
    Activity.Nocturnal==0 & Activity.Crepuscular==1 & Activity.Diurnal==0~ "C",
    Activity.Nocturnal==0 & Activity.Crepuscular==0 & Activity.Diurnal==1~ "D",
    Activity.Nocturnal==1 & Activity.Crepuscular==0 & Activity.Diurnal==0~ "N")) %>%
  # select only var of interest
  select(binomial, Habitat, Body.mass, ln.Mass.g,
         Main.diet, Activity, For.niche,
         insular_endemic) %>%
  filter(!is.na(Habitat))

summary(data_ft_mam %>% mutate_if(is.character, as.factor) %>% 
          select(-binomial))

# some NAs because of imputed values or db inaccuracy
# NA for Activity

na_activity <- pull(data_ft_mam %>% filter(is.na(Activity)) %>%
                      mutate(binomial = as.character(binomial)), binomial)

missing_activity <- lapply(df_mam_5t_imputed, function(x){
  x %>% filter(binomial %in% na_activity) %>% 
    select(binomial, contains("Activity"))
})
missing_activity_df <- missing_activity[[1]]
for(i in 2:ncol(missing_activity_df)){
  missing_activity_df[,i] <- 
    apply(do.call("rbind",lapply(missing_activity,"[[", i)), 2, sum)
}
new_activity <- data.frame(
  binomial = missing_activity_df$binomial,
  Activity = colnames(missing_activity_df)[
    apply(missing_activity_df,1,which.max)]) %>%
  mutate(Activity = case_when(
    Activity == "Activity.Diurnal" ~ "D",
    Activity == "Activity.Nocturnal" ~ "N"))

# NA for For.niche
na_niche <- pull(data_ft_mam %>% filter(is.na(For.niche)) %>%
                   mutate(binomial = as.character(binomial)), binomial)

missing_niche <- lapply(df_mam_5t_imputed, function(x){
  x %>% filter(binomial %in% na_niche) %>%
    select(binomial, `1:m_nbAr`, `1:m_nbG`, `1:m_nbM`,
           `1:m_nbA`, `1:m_nbS`) %>%
    rename(Arboreal = `1:m_nbAr`,
           Ground = `1:m_nbG`,
           Marine = `1:m_nbM`,
           Aerial = `1:m_nbA`,
           Scansorial = `1:m_nbS`)
})
missing_niche_df <- missing_niche[[1]]
for(i in 2:ncol(missing_niche_df)){
  missing_niche_df[,i] <-
    apply(do.call("rbind",lapply(missing_niche,"[[", i)), 2, sum)
}
# missing_niche_df

new_niche <- data.frame(
  binomial = missing_niche_df$binomial,
  For.niche = colnames(missing_niche_df)[
    apply(missing_niche_df,1,which.max)])

# take the modality with the highest sum 

for(i in 1:nrow(data_ft_mam)){
  for (j in 1:nrow(new_niche)){
    if (data_ft_mam$binomial[[i]]==new_niche$binomial[[j]]){
      data_ft_mam$For.niche[[i]] <- new_niche$For.niche[[j]]
    }
  }
  for (k in 1:nrow(new_activity)){
    if (data_ft_mam$binomial[[i]]==new_activity$binomial[[k]]){
      data_ft_mam$Activity[[i]] <- new_activity$Activity[[k]]
    }
  }
}
data_ft_mam <- data_ft_mam %>% 
  column_to_rownames(binomial) %>%
  mutate_if(is.character, as.factor)

summary(data_ft_mam %>% select(-binomial))

# -------------------------------- REPTILES 

df_out_rept <- df_rept_5t_imputed[[1]]
for(i in 2:ncol(df_out_rept)){
  df_out_rept[,i] <- 
    apply(do.call("rbind",lapply(df_rept_5t_imputed,"[[", i)), 2, median)
} 

colSums(is.na(df_out_rept))

data_ft_rept <- df_out_rept %>%
  # select only insular endemic species
  # filter(insular_endemic==1) %>% # KEEP ALL NATIVE SP 
  # Habitat trait
  mutate(hab_sum=as.character(hab_sum)) %>%
  mutate(hab_sum=if_else(hab_sum>4,"5",hab_sum)) %>%
  mutate(Habitat=factor(hab_sum, ordered = T,
                        levels = c("1","2","3","4","5"))) %>%
  # Body mass
  mutate(Body.mass=factor(ntile(max_log10_BM_g,5), ordered = T,
                          levels = c("1","2","3","4","5"))) %>%
  # Reproductive mode
  mutate(Repro.mode = case_when(
    `1:m_nbMixed`==1 ~ "Mixed",
    `1:m_nbOviparous`==1 ~ "Oviparous",
    `1:m_nbViviparous`==1 ~ "Viviparous")) %>%
  # Activity period
  mutate(Activity = case_when(
    `1:m_nbCathemeral`==1 ~ "Cathemeral",
    `1:m_nbDiurnal`==1 ~ "Diurnal",
    `1:m_nbNocturnal`==1 ~ "Nocturnal")) %>%
  # Substrate
  # consider unique substrate (Arb, Fos, Sax, Aqu, Ter)
  # multiple niches 1) involving aquatic niche (MAqu) or 2) not aqu (MNoAqu)
  mutate(For.niche = case_when(
    Sub_Fossorial==1 & Sub_Terrestrial==0 & Sub_Semi_Aquatic==0 & 
      Sub_Arboreal==0 & Sub_Saxicolous==0 ~ "Fos",
    Sub_Fossorial==0 & Sub_Terrestrial==1 & Sub_Semi_Aquatic==0 & 
      Sub_Arboreal==0 & Sub_Saxicolous==0 ~ "Ter",
    Sub_Fossorial==0 & Sub_Terrestrial==0 & Sub_Semi_Aquatic==1 & 
      Sub_Arboreal==0 & Sub_Saxicolous==0 ~ "Aqu",
    Sub_Fossorial==0 & Sub_Terrestrial==0 & Sub_Semi_Aquatic==0 & 
      Sub_Arboreal==1 & Sub_Saxicolous==0 ~ "Arb",
    Sub_Fossorial==0 & Sub_Terrestrial==0 & Sub_Semi_Aquatic==0 & 
      Sub_Arboreal==0 & Sub_Saxicolous==1 ~ "Sax",
    Sub_Arboreal+Sub_Fossorial+Sub_Saxicolous+Sub_Semi_Aquatic+
      Sub_Terrestrial>1 ~ "Multi")) %>%
  # select only var of interest
  select(binomial_iucn, Habitat, Body.mass, max_log10_BM_g,
         Repro.mode, Activity, For.niche, insular_endemic) %>%
  mutate_if(is.character, as.factor) %>%
  rename(binomial = binomial_iucn)

summary(data_ft_rept%>%select(-binomial))

# NA for Activity
na_activity <- pull(data_ft_rept %>% filter(is.na(Activity)) %>%
                      mutate(binomial = as.character(binomial)), binomial)

missing_activity <- lapply(df_rept_5t_imputed, function(x){
  x %>% filter(binomial_iucn %in% na_activity) %>% 
    select(binomial_iucn, `1:m_nbCathemeral`, 
           `1:m_nbDiurnal`, `1:m_nbNocturnal`) %>%
    rename(Cathemeral = `1:m_nbCathemeral`, 
           Diurnal = `1:m_nbDiurnal`, 
           Nocturnal = `1:m_nbNocturnal`)
})
missing_activity_df <- missing_activity[[1]]
for(i in 2:ncol(missing_activity_df)){
  missing_activity_df[,i] <- 
    apply(do.call("rbind",lapply(missing_activity,"[[", i)), 2, sum)
}

new_activity <- data.frame(
  binomial = missing_activity_df$binomial_iucn,
  Activity = colnames(missing_activity_df)[
    apply(missing_activity_df,1,which.max)])


# NA for For.niche
na_niche <- pull(data_ft_rept %>% filter(is.na(For.niche)) %>%
                   mutate(binomial = as.character(binomial)), binomial)

missing_niche <- lapply(df_rept_5t_imputed, function(x){
  x %>% filter(binomial_iucn %in% na_niche) %>%
    select(binomial_iucn, contains("Sub"))
})
missing_niche_df <- missing_niche[[1]]
for(i in 2:ncol(missing_niche_df)){
  missing_niche_df[,i] <-
    apply(do.call("rbind",lapply(missing_niche,"[[", i)), 2, sum)
}
# take the modality with the highest sum 
new_for_niche <- data.frame(
  binomial = missing_niche_df$binomial_iucn,
  For.niche = colnames(missing_niche_df)[
    apply(missing_niche_df,1,which.max)]) %>%
  mutate(For.niche = substr(sub(".*_", "", For.niche), 1, 3))

# NA for Repro.mode
na_repro <- pull(data_ft_rept %>% filter(is.na(Repro.mode)) %>%
                   mutate(binomial = as.character(binomial)), binomial)

missing_repro <- lapply(df_rept_5t_imputed, function(x){
  x %>% filter(binomial_iucn %in% na_repro) %>%
    select(binomial_iucn, `1:m_nbMixed`,
           `1:m_nbOviparous`, `1:m_nbViviparous`) %>%
    rename(Mixed = `1:m_nbMixed`,
           Oviparous = `1:m_nbOviparous`,
           Viviparous = `1:m_nbViviparous`)
})
missing_repro_df <- missing_repro[[1]]
for(i in 2:ncol(missing_repro_df)){
  missing_repro_df[,i] <-
    apply(do.call("rbind",lapply(missing_repro,"[[", i)), 2, sum)
}
# take the modality with the highest sum 
new_repro <- data.frame(
  binomial = missing_repro_df$binomial_iucn,
  Repro.mode = colnames(missing_repro_df)[
    apply(missing_repro_df,1,which.max)])

# Replace NAs in df_rept for all traits with missing values

for(i in 1:nrow(data_ft_rept)){
  for (j in 1:nrow(new_activity)){
    if (data_ft_rept$binomial[[i]]==new_activity$binomial[[j]]){
      data_ft_rept$Activity[[i]] <- new_activity$Activity[[j]]
    }}
  for (k in 1:nrow(new_for_niche)){
    if (data_ft_rept$binomial[[i]]==new_for_niche$binomial[[k]]){
      data_ft_rept$For.niche[[i]] <- new_for_niche$For.niche[[k]]
    }}
  for (l in 1:nrow(new_repro)){
    if (data_ft_rept$binomial[[i]]==new_repro$binomial[[l]]){
      data_ft_rept$Repro.mode[[i]] <- new_repro$Repro.mode[[l]]
    }}
}

summary(data_ft_rept %>% select(-binomial))
data_ft_rept <- data_ft_rept %>% 
  mutate(binomial = as.character(binomial))




# ---------     Create assemblage categories (threat)    -------------

# For each class, create assemblage
# bind df_all and threat to create IAS-T/NT var (and other-T/NT)
# select only insular tetrapods
df_all_threat <- left_join(df_all, all_threats, by="binomial")
colSums(is.na(df_all_threat))
df_all_threat[is.na(df_all_threat)] = 0

df_all_threat <- df_all_threat %>%
  mutate(threatened = if_else(category %in% c("CR","EN","VU"), 1, 0)) %>%
  mutate(ias_8.1 = if_else(category=="DD", 0, ias_8.1),
         other_threat = if_else(category=="DD", 0, other_threat)) %>%
  mutate(threatened_ias = paste(threatened, ias_8.1, other_threat, sep="_")) %>%
  mutate(group = case_when(
    threatened_ias == "1_1_0" ~ "IAS_T",
    threatened_ias == "0_1_0" ~ "IAS_NT",
    threatened_ias == "1_0_1" ~ "other_T",
    threatened_ias == "0_0_1" ~ "other_NT",
    threatened_ias == "1_0_0" | threatened_ias == "0_0_0"~ "not_included"
  )) 
table(df_all_threat$group)
# check for all IAS-T/NT (ias_8.1 = 1) to be in CR EN LC NT VU cate (not in DD)
table(df_all_threat$category, df_all_threat$group)


# ---------------   Save output files   -----------------------

setwd("Z:/THESE/6_Projects/predict_vulnerability/Output/Data_clean")

saveRDS(data_ft_mam, "02_ft_all_native_mammals")
saveRDS(data_ft_rept, "02_ft_all_native_reptiles")
saveRDS(df_birds_ft, "02_ft_all_native_birds")
saveRDS(df_all_threat, "02_threats_all_native_tetrapods")

