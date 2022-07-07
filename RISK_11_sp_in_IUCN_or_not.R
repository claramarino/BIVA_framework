# all sp with native range


rm(list=ls())

library(tidyverse)
library(rredlist)

# open species list 

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
fold <- list.dirs(main_dir)

all_sp_list <- vector(mode = "list", length = length(fold)-1)

for(i in 1:length(all_sp_list)){
  fil <- list.files(fold[i+1])
  sp_file <- fil[grepl("sp_list", fil)]
  all_sp_list[[i]] <- readRDS(paste0(fold[i+1], "/", sp_file))
}

all_sp <- bind_rows(all_sp_list)


gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")


sp_gbif <- left_join(all_sp, gbif_taxo_2, by = "scientificname")

sum(duplicated(sp_gbif$species.x))

sp_unique <- tolower(unique(sp_gbif$species.x))

# check for synonyms first
# then retrieve iucn information

# syno_ias <- data.frame()
# for (i in 1:length(sp_unique)){
#   obj <- rl_synonyms(sp_unique[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   syno_ias <- bind_rows(syno_ias, obj$result)
#   saveRDS(syno_ias,"Output/Synonyms/11_IUCN_synonyms_ias_from_gbif")
# }

syno_ias <- readRDS("Output/Synonyms/11_IUCN_synonyms_ias_from_gbif") %>% 
  distinct(accepted_name, synonym) %>%
  mutate_all(tolower)

# iucn_search_ias <- data.frame()
# for (i in 1:length(sp_unique)){
#   obj <- rl_search(sp_unique[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   if(is.null(obj$result)){
#     acc_name = syno_ias$accepted_name[syno_ias$synonym==sp_unique[i]]
#     if(length(acc_name)==0){ next }
#     obj_acc <- data.frame()
#     for (j in length(acc_name)){
#       obj2 <- rl_search(acc_name[j],
#                        key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#       obj2$result$ias_name = sp_unique[i]
#       obj_acc <- bind_rows(obj_acc, obj2$result)
#     }
#     iucn_search_ias <- bind_rows(iucn_search_ias, obj_acc)
#   } else {
#     obj$result$ias_name = sp_unique[i]
#     iucn_search_ias <- bind_rows(iucn_search_ias, obj$result)
#   }
#   saveRDS(iucn_search_ias,"Output/Native_exotic_range/11_IAS_sp_in_IUCN")
# }

iucn_search_ias <- readRDS("Output/Native_exotic_range/11_IAS_sp_in_IUCN") 
ias_in_iucn <- iucn_search_ias %>%
  filter(!is.na(taxonid)) %>% 
  mutate(ias_name = tolower(ias_name)) %>% 
  distinct(ias_name) %>% pull(ias_name)

ias_not_in_iucn <- sp_unique[!sp_unique %in% ias_in_iucn]



