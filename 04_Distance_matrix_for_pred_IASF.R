# Evaluation of status prediction for IAS-F species
# USing distance matrix based on gower distance between species
# As in mechabirds project

rm(list=ls())

library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(tidyr)
library(cluster)
library(tibble)
library(ggplot2)

# Load trait data
data_ft_mam <- readRDS("Output/Data_clean/02_ft_all_native_mammals")
data_ft_rept <- readRDS("Output/Data_clean/02_ft_all_native_reptiles")
data_ft_bird <- readRDS("Output/Data_clean/02_ft_all_native_birds")
# load threat data (global pool)
df_all_threat <- readRDS("Output/Data_clean/02_threats_all_native_tetrapods")
# load severe native (conservative pool)
natives_to_include <- readRDS("Output/Data_clean/01_native_names_to_model")

df_all_threat_conserv <- df_all_threat %>%
  filter(Class!="Amphibia") %>%
  mutate(Assoc_ias = ifelse(binomial %in% natives_to_include$all,1, 0),
         Assoc_ias_2more = ifelse(binomial %in% natives_to_include$all2more,1, 0),
         Assoc_ias_severe = ifelse(binomial %in% natives_to_include$severe,1, 0),
         Assoc_ias_severe2more = ifelse(binomial %in% natives_to_include$severe2more,1, 0))

table(df_all_threat_conserv$Assoc_ias_severe,
      df_all_threat_conserv$group)
table(df_all_threat_conserv$Assoc_ias_severe2more,
      df_all_threat_conserv$group)
table(df_all_threat_conserv$Assoc_ias_2more,
      df_all_threat_conserv$group)
table(df_all_threat_conserv$group, df_all_threat_conserv$Class)


# --------------- Create a list with the 3 classes -----------------

data_ft_bird$binomial = rownames(data_ft_bird)
# traits to keep for each class 
traits_mam <- c("Habitat","ln.Mass.g", "Main.diet","Activity",
                "For.niche", "insular_endemic")
traits_rept <- c("Habitat","max_log10_BM_g", "Repro.mode","Activity", 
                 "For.niche", "insular_endemic")

# select good traits for each class
data_ft_mam_ok <- data_ft_mam %>%
  select(binomial, all_of(traits_mam))
data_ft_rept_ok <- data_ft_rept %>%
  select(binomial, all_of(traits_rept))
# keep same traits as in mechabirds for birds

# create list
data_ft_list <- list(
  Bird = data_ft_bird,
  Mam = data_ft_mam_ok,
  Rept = data_ft_rept_ok
)

# shape data_ft for calculating distance matrix
data_ft_list <- lapply(data_ft_list, function(x){
  rownames(x) <- x$binomial
  x %>% select(-binomial) %>%
    mutate_if(is.character, as.factor)
})


####---------------- Set up species groups ------------------####

assoc_threat_list <- list(
  Bird = df_all_threat_conserv %>% filter(Class=="Aves"), 
  Mam = df_all_threat_conserv %>% filter(Class=="Mammalia"), 
  Rept = df_all_threat_conserv %>% filter(Class=="Reptilia"))

# ias_threatened sp (not conservative)
ias_t_all_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(group == "IAS_T")})
# threatened sp severely assoc to a named ias  
ias_t_severe_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(Assoc_ias_severe == 1 & group == "IAS_T")})
# all species severely associated to a named IAS
ias_a_severe_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(Assoc_ias_severe == 1)})

  
####-------------- Compute distance matrix -----------------####

# distance matrix based on traits 
dis.gower_list <- lapply(data_ft_list, daisy)

lapply(dis.gower_list, summary)
# globalement, les distances entre les oiseaux sont plus faibles que dans les autres classes
# comme dans soares, sélectionner un pool de traits qui permet de maximiser les distances ?

matrix_list <- lapply(dis.gower_list, as.matrix)
lapply(matrix_list, dim)

# Evaluate the prediction using all IAS-T (or IAS-A) sp 
# as a group with a known association with IAS
# select only those species in rows & in columns of the matrix
# <=> group to compare species to 

# SELECT HERE THE GROUP TO COMPARE (IAS-T, IAS-A, IAS-A severe, etc.)
comp_group_list <- lapply(ias_t_all_list, function(x){
  
})

no_dd <- unique(as.character(pull(ei_no_dd, Species1)))

df_dist <- as.data.frame(matrix)  %>%
  select(all_of(no_dd)) %>%
  rownames_to_column("no_DD_sp") %>%
  filter(no_DD_sp %in% no_dd) %>%
  column_to_rownames("no_DD_sp")

# create species to group correspondence
corresp <- ei_no_dd %>% distinct(Species1, group2)

#### Method 'Out of the bag' ####

# fill the  table for each of the 574 species (to remove in the predictors)

#------------- Initialize tables for metrics

# for min and mean distance
long_df_all_sp <- data.frame(
  no_DD_sp = character(0), group2 = character(0),
  mean_dist = numeric(0), min_dist = numeric(0))
# 10 closest neighbors
k_closest <- k_closest_sp <- data.frame(matrix(ncol = 10, nrow = 0)) %>%
  mutate_all(as.character)
colnames(k_closest) <- colnames(k_closest_sp) <- paste("closest_", 1:10, sep="")
# species in buffers 
d_buffer <- c("d1" = 0.21, 
              "d2" = mean(as.matrix(df_dist))) 
buffer_df_all <- data.frame(
  group2 = character(0), mean_buffer = numeric(0), n_b = numeric(0),
  prop_b = numeric(0), no_DD_sp = character(0),d_buffer=character(0))
#source function for calculating sp in buffer, mean dist & prop in buffer
source("R/find_sp_in_buffer.R")

#------------ Loop on all species

for (sp in no_dd){
  
  # test iteration
  #sp = no_dd[6]
  
  # select focus species & remove it from the predictors
  df_dist_out <- df_dist %>%
    select(-all_of(sp))%>%
    rownames_to_column("no_DD_sp") %>%
    filter(no_DD_sp == sp) %>%
    column_to_rownames("no_DD_sp")
  
  #___________________ Min & mean distance with groups
  
  long_df <- df_dist_out %>% rownames_to_column("no_DD_sp") %>%
    pivot_longer(!no_DD_sp, names_to = "sp2", values_to = "dist")
  
  long_df_sp <- left_join(long_df, corresp %>% rename(sp2 = Species1), 
                          by = "sp2") %>%
    group_by(no_DD_sp, group2) %>%
    summarise(mean_dist = mean(dist), min_dist = min(dist))
  
  # implement in final long df with all no_dd species
  long_df_all_sp <- bind_rows(long_df_all_sp, long_df_sp)
  
  #____________________ 10 closest neighbours
  
  col1 <- as.data.frame(t(df_dist_out)) %>% rownames_to_column("sp")
  sorted <- col1[order(col1[,2]),]
  for (k in 1:10){
    k_closest_sp[1,k]<- sorted$sp[k]
  } 
  row.names(k_closest_sp) <- sp
  
  # implement in final k_closest dataset
  k_closest <- bind_rows(k_closest, k_closest_sp)
  
  #____________________ species in buffers
  
  buffer_d1 <- find_sp_in_buffer_gp(df_dist_out, sp, d_buffer[1])
  buffer_d2 <- find_sp_in_buffer_gp(df_dist_out, sp, d_buffer[2])
  
  buffer_d1[[1]] <- buffer_d1[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d1")
  buffer_d2[[1]] <- buffer_d2[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d2")
  
  buffer_df_all <- bind_rows(buffer_df_all, 
                             do.call(rbind.data.frame, buffer_d1),
                             do.call(rbind.data.frame, buffer_d2))
  
}

#___________________ Final metrics dataset

k_closest_lg <- k_closest %>%
  rownames_to_column("no_DD_sp") %>%
  pivot_longer(!no_DD_sp, names_to = "rank", values_to = "sp_2")

k_closest_gps <- left_join(k_closest_lg, corresp %>%
                             rename(sp_2 = Species1), 
                           by = "sp_2") %>% 
  group_by(no_DD_sp, group2) %>% summarise(n_10_closest = n())

buffer_d1_df <- buffer_df_all %>%
  filter(d_buffer == "d1") %>%
  rename(mean_b1 = mean_buffer,
         n_b1 = n_b,
         prop_b1 = prop_b) %>%
  select(-d_buffer)

buffer_d2_df <- buffer_df_all %>%
  filter(d_buffer == "d2") %>%
  rename(mean_b2 = mean_buffer,
         n_b2 = n_b,
         prop_b2 = prop_b) %>%
  select(-d_buffer)

evaluate_no_dd <- full_join(
  # Add k closest neighbours
  full_join(long_df_all_sp, k_closest_gps, by=c("no_DD_sp","group2")),
  # Add buffer metrics
  # join with d2 first because no NA
  full_join(buffer_d2_df, buffer_d1_df, by=c("no_DD_sp","group2"))) %>%
  # replace NA by 0 for proportions (but nor for means)
  mutate_at(c("n_10_closest","n_b1","prop_b1", "n_b2","prop_b2"), 
            ~ ifelse(is.na(.), 0, .))

#saveRDS(evaluate_no_dd, "Output/Eval_DD/21_metrics_evaluate_no_dd_ft_pc1to4")
saveRDS(evaluate_no_dd, "Output/Eval_DD/21_metrics_evaluate_no_dd_ft")


