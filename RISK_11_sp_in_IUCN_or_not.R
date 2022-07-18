# all sp with native range


rm(list=ls())

library(tidyverse)
library(rredlist)

# open species list 

all_ias_with_occ <- readRDS("Output/Occurrences_clean_taxo_ok/RISK_03_ias_list_with_occ")

sp_unique <- tolower(unique(all_ias_with_occ$new_species))

# check for synonyms first
# then retrieve iucn information

# syno_ias <- data.frame()
# for (i in 1:length(sp_unique)){
#   obj <- rl_synonyms(sp_unique[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   syno_ias <- bind_rows(syno_ias, obj$result)
#   saveRDS(syno_ias,"Output/Synonyms/RISK_11_IUCN_synonyms_ias_from_gbif")
# }

syno_ias <- readRDS("Output/Synonyms/RISK_11_IUCN_synonyms_ias_from_gbif") %>% 
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
#   saveRDS(iucn_search_ias,"Output/Native_exotic_range/RISK_11_IAS_sp_in_IUCN")
# }

iucn_search_ias <- readRDS("Output/Native_exotic_range/RISK_11_IAS_sp_in_IUCN") 

ias_in_iucn <- iucn_search_ias %>%
  filter(!is.na(taxonid)) %>% 
  mutate(ias_name = tolower(ias_name)) %>% 
  select(ias_name, scientific_name, kingdom:class, category, 
           population_trend:terrestrial_system) %>%
  distinct()


length(unique(iucn_search_ias$scientific_name))

table(ias_in_iucn$terrestrial_system)

table(ias_in_iucn$class)
table(ias_in_iucn$kingdom)
table(ias_in_iucn$phylum)


# cb de sp par classe à la base ?
gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")
sp_class <- left_join(all_ias_with_occ,
                      gbif_taxo_2 %>% mutate(new_species = species), 
                      by="new_species") %>%
  select(new_species, kingdom, phylum, class) %>% distinct()


table(sp_class$class)
table(sp_class$kingdom)
table(sp_class$phylum)



#### Save species list inIUCN/ not in IUCN for script 13 #####

ias_not_in_iucn <- sp_unique[!sp_unique %in% ias_in_iucn$ias_name]


sp_not_in_iucn_table <- all_ias_with_occ %>%
  mutate(spe_lower = tolower(new_species)) %>%
  filter(spe_lower %in% ias_not_in_iucn)

saveRDS(sp_not_in_iucn_table, "Output/Native_exotic_range/RISK_11_IAS_sp_NOT_in_IUCN")




all_ias_with_occ_iucn_check <- all_ias_with_occ%>%
  mutate(spe_lower = tolower(new_species)) %>%
  mutate(ias_in_iucn = if_else(spe_lower %in% ias_in_iucn$ias_name, "YES","NO"))

saveRDS(all_ias_with_occ_iucn_check, "Output/Native_exotic_range/RISK_11_ias_list_with_occ_IUCN_check")

