# all sp with native range


rm(list=ls())

library(tidyverse)
library(rredlist)

# open species list 

main_dir <- "Z:/THESE/6_Projects/predict_vulnerability/Output/Occurrences_clean/"
fold <- list.dirs(main_dir)
clean_fold <- fold[-c(1,5)]

clean_fold


all_sp_list <- vector(mode = "list", length = length(clean_fold))

for(i in 1:length(all_sp_list)){
  fil <- list.files(clean_fold[i])
  sp_file <- fil[grepl("sp_list", fil)]
  all_sp_list[[i]] <- readRDS(paste0(clean_fold[i], "/", sp_file))
}

all_sp <- bind_rows(all_sp_list)

length(unique(all_sp$specieskey))


# after cleaning, some species have less than ten occurrences
all_sp %>% filter(n_occ<10)

# due to taxo pbm (forget to look for some species but only search for synonyms?)

gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")

uk <- unique(gbif_taxo_2$usagekey)
spk <- unique(gbif_taxo_2$specieskey)

spk_in_uk <- spk[spk %in% uk]

spk_not_uk <- spk[!spk %in% uk]

tax_prob <- gbif_taxo_2 %>% filter(specieskey %in% spk_not_uk)



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


#### Separate species in IUCN for retrieve native range #####







#### Save species list not in IUCN for script 13 #####

ias_not_in_iucn <- sp_unique[!sp_unique %in% ias_in_iucn]

sp_not_in_iucn_table <- gbif_taxo_2 %>%
  mutate(spe_lower = tolower(species)) %>%
  filter(spe_lower %in% ias_not_in_iucn)


saveRDS(sp_not_in_iucn_table, "Output/Native_exotic_range/11_IAS_sp_NOT_in_IUCN")


