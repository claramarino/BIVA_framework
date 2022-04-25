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
    select(ias_8.1, IAS_T, no_IAS_T, Assoc_ias:Assoc_ias_2more_t, global)
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


# all species associated to a named IAS
ias_as_2more_list <- lapply(assoc_threat_list, function(x){
  x %>% filter(Assoc_ias_2more == 1) %>%
    pull(binomial)})







################################################################
#              Distance matrix for predicting IAS-as
################################################################

# test sur la matrice de distances entre les sp
lapply(dis.gower_list, summary)
# globalement, les distances entre les oiseaux sont plus faibles que dans les autres classes
# comme dans soares, sélectionner un pool de traits qui permet de maximiser les distances ?
matrix_list <- lapply(dis.gower_list, as.matrix)
# for birds, keep only species within groups
matrix_list$Bird <- matrix_list$Bird[assoc_threat_list$Bird$binomial,
                                     assoc_threat_list$Bird$binomial]

lapply(matrix_list, dim)
lapply(assoc_threat_list, dim)


hist(matrix_list$Bird, breaks = 100)
hist(matrix_list$Mam, breaks = 100)
hist(matrix_list$Rept, breaks = 100)


library(ks)
dens_b <- density(matrix_list$Bird)
dens_m <- density(matrix_list$Mam)
dens_r <- density(matrix_list$Rept)

plot(dens_b)
plot(dens_m)
plot(dens_r)

library(pastecs)



# plot densities of distances and turningpoints of densities distribution of dist

# Birds
loc.min.max_b <- turnpoints(dens_b$y)
dens_b$x[loc.min.max_b$pits]
dens_b$y[loc.min.max_b$pits]
hist(matrix_list$Bird, breaks = 100)
plot(dens_b)
points(dens_b$x[loc.min.max_b$pits],
       dens_b$y[loc.min.max_b$pits],col="green")
abline(v = dens_b$x[loc.min.max_b$pits], col = "gray", lty = 3)

# Mammals
loc.min.max_m <- turnpoints(dens_m$y)
dens_m$x[loc.min.max_m$pits]
dens_m$y[loc.min.max_m$pits]
hist(matrix_list$Mam, breaks = 100)
plot(dens_m)
points(dens_m$x[loc.min.max_m$pits],
       dens_m$y[loc.min.max_m$pits],col="green")
abline(v = dens_m$x[loc.min.max_m$pits], col = "gray", lty = 3)

# Reptiles
loc.min.max_r <- turnpoints(dens_r$y)
dens_r$x[loc.min.max_r$pits]
dens_r$y[loc.min.max_r$pits]
hist(matrix_list$Rept, breaks = 100)
plot(dens_r)
points(dens_r$x[loc.min.max_r$pits],
       dens_r$y[loc.min.max_r$pits],col="green")
abline(v = dens_r$x[loc.min.max_r$pits], col = "gray", lty = 3)


dist_buffer <- list(
  Bird = 0.097282588, 
  Mam = 0.15126421,
  Rept = 0.15231439)



# Evaluate the prediction using all IAS-T (or IAS-A) sp 
# as a group with a known association with IAS
# select only those species in rows & in columns of the matrix
# <=> group to compare species to 


# Initialize output df
long_df_all_sp_bmr <- list(
  Bird = data.frame(),
  Mam = data.frame(),
  Rept = data.frame())

for(class in 1:3){
  
  # select ias_as species for class i
  ias_as <- ias_as_2more_list[[class]]
  
  # select group to compare
  comp_group <- assoc_threat_list[[class]] %>% 
    select(binomial, Assoc_ias_2more) %>%
    rename(sp2 = binomial)
  
  df_dist <- as.data.frame(matrix_list[[class]])%>%
    rownames_to_column("ias_as_sp") %>%
    filter(ias_as_sp  %in% ias_as) %>%
    column_to_rownames("ias_as_sp")
  
  for (sp in ias_as){
    sp = ias_as[15]
    
    # filter distances only for those species 
    df_dist_out <- df_dist  %>%
      select(-all_of(sp))%>%
      rownames_to_column("ias_as_sp") %>%
      filter(ias_as_sp == sp) %>%
      column_to_rownames("ias_as_sp")
    
    long_df <- df_dist_out %>% rownames_to_column("ias_as_sp") %>%
      pivot_longer(!ias_as_sp, names_to = "sp2", values_to = "dist")
    
    long_df_sp <- left_join(long_df, comp_group, 
                            by = "sp2") %>%
      group_by(ias_as_sp, Assoc_ias_2more) %>%
      summarise(mean_dist = mean(dist), min_dist = min(dist))
    
    # implement in final long df with all ias-as species
    long_df_all_sp_bmr[[class]] <- bind_rows(long_df_all_sp_bmr[[class]], 
                                             long_df_sp)
    
    
    
    # test density
    # dens <- density(long_df$dist)
    # 
    # loc.min.max <- turnpoints(dens$y)
    # bw <- dens$bw
    # md <- dens$x[loc.min.max$pits][1]
    # 
    # long_df_crit <- long_df %>%
    #   mutate(band_width = bw - dist, 
    #          min_tp = md - dist)
    # 
    # ggplot(left_join(long_df, comp_group, by = "sp2") %>%
    #          mutate(group = as.character(Assoc_ias_2more)), aes(x = dist)) +
    #   geom_density(aes(group = group, fill = group), alpha = 0.5)
    }
  
  
}
saveRDS(long_df_all_sp_bmr, "Output/distance_traits/04_dist_ias_as_to_group_eval")


eval_pred <- lapply(long_df_all_sp_bmr, function(x){
  
  status_all_meth <- data.frame(matrix(ncol = 3, nrow = length(unique(x$ias_as_sp))))
  names(status_all_meth) <- c("ias_as_sp","min_dist","mean_dist")
  status_all_meth$ias_as_sp <- unique(x$ias_as_sp)
  
  vars_min <- c("min_dist","mean_dist")
  
  for(sp in status_all_meth$ias_as_sp){
    
    evaluate_sp <- x %>% filter(ias_as_sp == sp)
    
    for(col in vars_min){
      status_all_meth[which(status_all_meth$ias_as_sp==sp),col] <-
        evaluate_sp$Assoc_ias_2more[which.min(pull(evaluate_sp,col))]}
  }
  
  return(status_all_meth)
  
})

lapply(eval_pred, function(x){
  colSums(x %>% select(min_dist, mean_dist))/nrow(x)})
















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


