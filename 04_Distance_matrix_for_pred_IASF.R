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

# Load distance matrix

dis.gower_list <- readRDS("Output/03_distance_matrix_bmr")

# Load assemblages
df_all_threat_conserv <- readRDS("Output/Data_clean/03_df_all_threat_assoc_bmr")
assoc_threat_list <- readRDS("Output/Data_clean/03_df_all_threat_assoc_bmr_with_traits")




assembl_threat_list <- lapply(assoc_threat_list, function(x){
  rownames(x) = x$binomial
  to_transpose  <- x %>% 
    mutate(IAS_T = if_else(group=="IAS_T", 1, 0),
           global = 1,
           no_IAS_T = global - IAS_T) %>%
    select(ias_8.1, IAS_T, no_IAS_T, Assoc_ias:Assoc_ias_severe2more, global)
  transposed <- t(as.matrix(to_transpose))
  return(transposed)
})

# ias_threatened sp (not conservative)
ias_t_all_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(group == "IAS_T")})
# threatened sp severely assoc to a named ias  
ias_t_severe_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(Assoc_ias_severe == 1 & group == "IAS_T")})
# all species severely associated to a named IAS
ias_a_severe_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(Assoc_ias_severe == 1)})

################################################################
#              Distance matrix for predicting IAS-as
################################################################

# test sur la matrice de distances entre les sp

dis.gower_list <- lapply(data_ft_list, cluster::daisy)
lapply(dis.gower_list, summary)
# globalement, les distances entre les oiseaux sont plus faibles que dans les autres classes
# comme dans soares, sélectionner un pool de traits qui permet de maximiser les distances ?
matrix_list <- lapply(dis.gower_list, as.matrix)
lapply(matrix_list, dim)

hist(matrix_list$Bird, breaks = 100)
hist(matrix_list$Mam, breaks = 100)
hist(matrix_list$Rept, breaks = 100)


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


