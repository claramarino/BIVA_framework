rm(list=ls())


# Load packages
library(tidyverse)
library(BAT)
library(alphahull)
library(ape)
library(hypervolume)
# # Automatically install required packages, which are not yet in library
# packages <- c("tidyverse","BAT")
# new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages); rm(new.packages)
# # Load required packages
# l <- sapply(packages, require, character.only = TRUE); rm(packages, l)


# Load trait data
data_ft_mam <- readRDS("Output/Data_clean/02_ft_all_native_mammals")
data_ft_rept <- readRDS("Output/Data_clean/02_ft_all_native_reptiles")
data_ft_bird <- readRDS("Output/Data_clean/02_ft_all_native_birds")
# load threat data (global pool)
df_all_threat <- readRDS("Output/Data_clean/02_threats_all_native_tetrapods")
# load severe native (conservative pool)
natives_to_include <- readRDS("Output/Data_clean/01_native_names_to_model")

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
rm(data_ft_bird, data_ft_mam, data_ft_mam_ok, data_ft_rept, data_ft_rept_ok)

# shape data_ft for calculating distance matrix
data_ft_list <- lapply(data_ft_list, function(x){
  rownames(x) <- x$binomial
  x %>% select(-binomial) %>%
    mutate_if(is.character, as.factor)
})


# --------------- Functional space axes -----------------

# Build functional space axis with pcoa 
# Same method as in mechabirds
# distance matrix based on traits 
# 
# dis.gower_list <- lapply(data_ft_list, cluster::daisy)
# saveRDS(dis.gower_list, "Output/03_distance_matrix_bmr")
# pcoa_list <- lapply(dis.gower_list, ape::pcoa) # very long (30 min +)
# # save in RDS
# saveRDS(pcoa_list, "Output/03_pcoa_coordinates_all_tetrapods")
pcoa_list <- readRDS("Output/03_pcoa_coordinates_all_tetrapods")

# keep best number of axes
nbdim = 15
mat_coord_list <- lapply(pcoa_list, function(x){
  nbdim_select<-min(nbdim,ncol(x$vectors) )
  # keeping species coordinates on the 'nbdim' axes
  mat_coord<-x$vectors[,1:nbdim_select]
  return(mat_coord)
  })

mat_coord_best <- mat_coord_list

for(i in names(data_ft_list)){
  row.names(mat_coord_list[[i]])<-row.names(data_ft_list[[i]])
  colnames(mat_coord_list[[i]])<-paste("PC",1:ncol(mat_coord_list[[i]]),sep="")
  print(i)
  
  # check number of axis to retain 80% of variance
  sum_var = 0
  k=0
  while (sum_var < 0.8){
    k = k+1
    sum_var = sum_var + pcoa_list[[i]]$values$Eigenvalues[k]/
      sum(pcoa_list[[i]]$values$Eigenvalues[which(pcoa_list[[i]]$values$Eigenvalues>0)])
  print(sum_var)
  }
  mat_coord_best[[i]] <- as.data.frame(mat_coord_list[[i]][,1:k])
}


#----------- Assemblage IAS-T, IAS-A, IAS --------------------


df_all_threat_conserv <- df_all_threat %>%
  filter(Class!="Amphibia") %>%
  mutate(Assoc_ias = ifelse(binomial %in% natives_to_include$all,1, 0),
         Assoc_ias_2more = ifelse(binomial %in% natives_to_include$all2more,1, 0),
         Assoc_ias_severe = ifelse(binomial %in% natives_to_include$severe,1, 0),
         Assoc_ias_severe2more = ifelse(binomial %in% natives_to_include$severe2more,1, 0))
#saveRDS(df_all_threat_conserv, "Output/Data_clean/03_df_all_threat_assoc_bmr")


assoc_threat_list <- list(
  Bird = df_all_threat_conserv %>% filter(binomial %in% rownames(data_ft_list$Bird)), 
  Mam = df_all_threat_conserv %>% filter(binomial %in% rownames(data_ft_list$Mam)), 
  Rept = df_all_threat_conserv %>% filter(binomial %in% rownames(data_ft_list$Rept)))

#saveRDS(assoc_threat_list, "Output/Data_clean/03_df_all_threat_assoc_bmr_with_traits")

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


# Objectif de la méthode :
# créer un hypervolume avec toutes les espèces IAS-associated (2 more)
# sélectionner les zones coeurs de l'hv (hotspots)
# évaluer quels points des espèces non associées tombent dans les zones coeur de l'hv
# répéter la procédure pour plusieurs hv (100) car probalistic doinc variable
# ne garder que les espèces qui tombent à chaque fois dans la zone coeur

#####################################################################
#              Evaluation of the method
#####################################################################

# Pour évaluer la méthode :
# tester si les IAS-associated tombent à chaque fois dans l'hv
# Pour chaque IAS-as, construire un hv sans elle (mais avec tout le reste du groupe)
# évaluer si le point tombe dans l'hv


pred_in_core <- vector(mode = "list", length = 3)
names(pred_in_core) <- names(assembl_threat_list)

# very long : from 5hours for mammals to 15h for half of birds
# do it in parallel for each group and separate birds in two

# for (i in 1:3){
#   print(names(assembl_threat_list)[i])
#   
#   # select assemblage of class i
#   assembl <- assembl_threat_list[[i]]
#   # select coordinates of all species that are in assemblage of class i
#   mat_coord_ok <- mat_coord_best[[i]][
#     pull(assoc_threat_list[[i]], binomial),]
#   
#   # set up species list 
#   sp_list <- pull(assoc_threat_list[[i]] %>% filter(Assoc_ias_2more == 1), 
#                   binomial)
#   
#   #### For birds separate in two list for parallelizing
#   # sp_list_b1 <- sp_list[1:(length(sp_list)/2)]
#   # sp_list_b2 <- sp_list[((length(sp_list)/2)+1):length(sp_list)]
#   
#   pred_in_core[[i]] <- data.frame(
#     binomial = sp_list,
#     is_in_core_0.01 = character(length(sp_list)),
#     is_in_core_0.05 = character(length(sp_list)))
#   
#   for (sp in sp_list){
#     print(paste0(which(sp_list == sp), "/", length(sp_list)))
#     
#     group_no_sp <- assembl["Assoc_ias_2more", ]
#     group_no_sp[sp] <- 0
#     
#     hv_without_sp <- kernel.build(group_no_sp, 
#                                   mat_coord_ok)
#     
#     point_to_test <- mat_coord_ok[sp,] 
#   
#     # for p = 0.01
#     hotspot_0.01 <- kernel.hotspots(hv_without_sp, prop = 0.01)
#     points_in_hotspot_0.01 <- hypervolume_inclusion_test(
#       hotspot_0.01, point_to_test, fast.or.accurate = "accurate")
# 
#     # for p = 0.05
#     hotspot_0.05 <- kernel.hotspots(hv_without_sp, prop = 0.05)
#     points_in_hotspot_0.05 <- hypervolume_inclusion_test(
#       hotspot_0.05, point_to_test, fast.or.accurate = "accurate")
#     
#     # store results
#     
#     pred_in_core[[i]]$is_in_core_0.01[pred_in_core[[i]]$binomial==sp] <- 
#       as.character(points_in_hotspot_0.01)
#     pred_in_core[[i]]$is_in_core_0.05[pred_in_core[[i]]$binomial==sp] <- 
#       as.character(points_in_hotspot_0.05)
#     
#   }
#   
#   
# }
# 
# saveRDS(pred_in_core$Mam, "Output/hv/03_ias_as_sp_in_core_mam")
# saveRDS(pred_in_core$Rept, "Output/hv/03_ias_as_sp_in_core_rept")
# saveRDS(pred_in_core$Bird, "Output/hv/03_ias_as_sp_in_core_bird_b1")
# saveRDS(pred_in_core$Bird, "Output/hv/03_ias_as_sp_in_core_bird_b2")
# 

#-------------------- Read files

pred_in_core <- list(
  Bird = bind_rows(
    readRDS("Output/hv/03_ias_as_sp_in_core_bird_b1"),
    readRDS("Output/hv/03_ias_as_sp_in_core_bird_b2")),
  Mam = readRDS("Output/hv/03_ias_as_sp_in_core_mam"),
  Rept = readRDS("Output/hv/03_ias_as_sp_in_core_rept")
)

lapply(pred_in_core, function(x){
  sum(x$is_in_core_0.01 == "TRUE")
  sum(x$is_in_core_0.05 == "TRUE")
})







#####################################################################
#              Find species in ias hv core
#####################################################################

#test sur 100 hv 

#------------- For mammals --------------

sp_in_hotspot_as_0.05 <- vector(mode = "list", length = 100)
sp_in_hotspot_as_0.01 <- vector(mode = "list", length = 100)

for (i in 1:100){
  #compute hv
  hv_ias_as_mam <- kernel.build(assembl_threat_list$Mam["Assoc_ias_2more",], 
                                mat_coord_best$Mam)
  # for p = 0.01
  hotspot_ias_as_0.01 <- kernel.hotspots(hv_ias_as_mam, prop = 0.01)
  points_in_hotspot_0.01 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.01, points_to_test_as,fast.or.accurate = "accurate")
    sp_in_hotspot_as_0.01[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.01,])
  
  # for p = 0.05
  hotspot_ias_as_0.05 <- kernel.hotspots(hv_ias_as_mam, prop = 0.05)
  points_in_hotspot_0.05 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.05, points_to_test_as,fast.or.accurate = "accurate")
  sp_in_hotspot_as_0.05[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.05,])
  
  saveRDS(sp_in_hotspot_as_0.01, "Output/hv/03_sp_in_hotspot_a2_0.01_mam")
  saveRDS(sp_in_hotspot_as_0.05, "Output/hv/03_sp_in_hotspot_a2_0.05_mam")
  
  print(i)
}





#--------------- for birds --------------

# test sur 100 hv 

sp_in_hotspot_as_0.05 <- vector(mode = "list", length = 100)
sp_in_hotspot_as_0.01 <- vector(mode = "list", length = 100)
# select sp that are not IAS-AS
non_ias_as_bird <- pull(assoc_threat_list$Bird %>% filter(Assoc_ias_severe != 1), 
                        binomial)
points_to_test_as <- mat_coord_best$Bird[non_ias_as_bird,] 


for (i in 1:100){
  #compute hv
  hv_ias_as_bird <- kernel.build(assembl_threat_list$Bird["Assoc_ias_severe",], 
                                 mat_coord_best$Bird[
                                   pull(assoc_threat_list$Bird, binomial),])
  # for p = 0.01
  hotspot_ias_as_0.01 <- kernel.hotspots(hv_ias_as_bird, prop = 0.01)
  points_in_hotspot_0.01 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.01, points_to_test_as,fast.or.accurate = "accurate")
  sp_in_hotspot_as_0.01[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.01,])
  
  # for p = 0.05
  hotspot_ias_as_0.05 <- kernel.hotspots(hv_ias_as_bird, prop = 0.05)
  points_in_hotspot_0.05 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.05, points_to_test_as,fast.or.accurate = "accurate")
  sp_in_hotspot_as_0.05[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.05,])
  
  print(i)
  
  saveRDS(sp_in_hotspot_as_0.01, "Output/hv/03_sp_in_hotspot_as_0.01_bird")
  saveRDS(sp_in_hotspot_as_0.05, "Output/hv/03_sp_in_hotspot_as_0.05_bird")
}

#--------------- for reptiles --------------

# test sur 100 hv 

sp_in_hotspot_as_0.05 <- vector(mode = "list", length = 100)
sp_in_hotspot_as_0.01 <- vector(mode = "list", length = 100)
# select sp that are not IAS-AS
non_ias_as_rept <- pull(assoc_threat_list$Rept %>% filter(Assoc_ias_severe != 1), 
                        binomial)
points_to_test_as <- mat_coord_best$Rept[non_ias_as_rept,] 


for (i in 1:100){
  #compute hv
  hv_ias_as_rept <- kernel.build(assembl_threat_list$Rept["Assoc_ias_severe",], 
                                 mat_coord_best$Rept)
  # for p = 0.01
  hotspot_ias_as_0.01 <- kernel.hotspots(hv_ias_as_rept, prop = 0.01)
  points_in_hotspot_0.01 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.01, points_to_test_as,fast.or.accurate = "accurate")
  sp_in_hotspot_as_0.01[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.01,])
  
  # for p = 0.05
  hotspot_ias_as_0.05 <- kernel.hotspots(hv_ias_as_rept, prop = 0.05)
  points_in_hotspot_0.05 <- hypervolume_inclusion_test(
    hotspot_ias_as_0.05, points_to_test_as,fast.or.accurate = "accurate")
  sp_in_hotspot_as_0.05[[i]] <- rownames(points_to_test_as[points_in_hotspot_0.05,])
  
  saveRDS(sp_in_hotspot_as_0.01, "Output/hv/03_sp_in_hotspot_as_0.01_rept")
  saveRDS(sp_in_hotspot_as_0.05, "Output/hv/03_sp_in_hotspot_as_0.05_rept")
  
  print(i)
}





#############################################################
#                         TESTS
#############################################################


# build hypervolume
# from BAT package (Soares)

hv_ias_t_mam <- kernel.build(assembl_threat_list$Mam[2,], mat_coord_best$Mam)
plot(hv_ias_t_mam)

rownames(assembl_threat_list$Mam)

hv_ias_as_mam <- kernel.build(assembl_threat_list$Mam["Assoc_ias_severe",], 
                              mat_coord_best$Mam)

# select sp that are not IAS-T or IAS-AS
non_ias_t_mam <- pull(assoc_threat_list$Mam %>% filter(group !="IAS_T"), binomial)
non_ias_as_mam <- pull(assoc_threat_list$Mam %>% filter(Assoc_ias_severe != 1), binomial)


points_to_test <- mat_coord_best$Mam[non_ias_t_mam,]
points_to_test_as <- mat_coord_best$Mam[non_ias_as_mam,] 

points_in_hv <- hypervolume_inclusion_test(hv_ias_t_mam, points_to_test,
                                           fast.or.accurate = "accurate")

points_falling_in <- points_to_test[points_in_hv,]
sum(points_in_hv)
length(points_in_hv)

hv_no_iast_mam <- kernel.build(assembl_threat_list$Mam["no_IAS_T",], mat_coord_best$Mam)

saveRDS(hv_no_iast_mam, "Output/hv/03_hv_no_iast_mam")
saveRDS(hv_ias_t_mam, "Output/hv/03_hv_iast_mam")


intersect_iast_no_iast_mam <- hypervolume_set(hv_ias_t_mam, hv_no_iast_mam,
                                              check.memory = F)

summary(intersect_iast_no_iast_mam)
get_volume(intersect_iast_no_iast_mam)
kernel.alpha(intersect_iast_no_iast_mam)


hotspot_ias_t <- kernel.hotspots(hv_ias_t_mam, prop = 0.05)
points_in_hotspot <- hypervolume_inclusion_test(hotspot_ias_t, points_to_test,
                                                fast.or.accurate = "accurate")
sum(points_in_hotspot)

data_intersect <- intersect_iast_no_iast_mam@HVList$Intersection@Data
# no data

hotspot_ias_as <- kernel.hotspots(hv_ias_as_mam, prop = 0.01)
points_in_hotspot_as <- hypervolume_inclusion_test(hotspot_ias_as, points_to_test_as,
                                                   fast.or.accurate = "accurate")
sum(points_in_hotspot_as)

# Test functions from mFD

library(mFD)

fric <- alpha.fd.multidim(
  sp_faxes_coord = as.matrix(mat_coord_best$Bird[pull(assoc_threat_list$Bird, binomial),]),
  asb_sp_w = assembl_threat_list$Bird,
  ind_vect = "fric")

plots_alpha <- alpha.multidim.plot(
  output_alpha_fd_multidim = fric,
  ind_nm = "fric",
  plot_asb_nm = "Assoc_ias_2more")

plots_alpha$"fric"$"patchwork"
