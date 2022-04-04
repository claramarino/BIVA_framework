# Focus on species associated to native amniotes
# Count IAS associated to each group/all groups

rm(list=ls())

library(tidyr)
library(dplyr)
library(rredlist)
library(tidyverse)


ias_x_native <- readRDS("Output/00_IAS_x_native_species")

# remove native amphibians
# remove unspecified sp

ias_x_native_spe <- ias_x_native %>%
  filter(Class!="Amphibia") %>%
  filter(!grepl("unspecified",ias_lower))

# number of IAS in total (on rept, mam, birds)
length(unique(ias_x_native_spe$ias_lower))

# number of IAS that are causing declines
severe_inter <- c("Very Rapid Declines", "Slow, Significant Declines", 
                  "Rapid Declines", "Causing/Could cause fluctuations")

length(unique(ias_x_native_spe %>%
         filter(severity %in% severe_inter) %>%
         pull(ias_lower)))

# number of IAS that have more than 1 interaction


# sp IAS that threatened 2-3 classes
ias_class <- as.data.frame(
  table(ias_x_native_spe$ias_lower, ias_x_native_spe$Class)) %>%
  pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
  rename(ias_name = Var1) %>%
  mutate(Total_inter = Aves + Mammalia + Reptilia) %>%
  mutate(Specificity = case_when(
    Aves == 0 & Mammalia == 0 ~ "Only_rept",
    Reptilia == 0 & Mammalia == 0 ~ "Only_Aves",
    Aves == 0 & Reptilia == 0 ~ "Only_mam",
    Aves > 0 & Mammalia > 0 ~ "Aves_mam",
    Reptilia > 0 & Mammalia > 0 ~ "Mam_rept",
    Aves > 0 & Reptilia > 0 ~ "Aves_rept")) %>%
  mutate(Specificity = if_else(Aves > 0 & Reptilia > 0 & Mammalia > 0,
                               "All_classes", Specificity))

table(ias_class$Specificity)
table(ias_class$Specificity[ias_class$Total_inter == 1])

length( ias_class$Total_inter[ias_class$Total_inter > 1] )
sum( ias_class$Total_inter[ias_class$Total_inter > 1] )

ias_class %>%
  group_by(Specificity) %>%
  summarise(nb_inter = sum(Total_inter))

names_2_or_more <- as.character(ias_class$ias_name[ias_class$Total_inter > 1])
# Find synonyms using rredlist and APi key to have all synonyms from IUCN 
# some species might be duplicated because of their synonyms
# ex sus domesticus = sus scrofa ?

# syno_ias <- data.frame()
# for (i in 1:length(names_2_or_more)){
#   obj <- rl_synonyms(names_2_or_more[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   syno_ias <- bind_rows(syno_ias, obj$result)
#   saveRDS(syno_ias,"Output/Synonyms/01_synonyms_ias_2_or_more")
# }

syno_ias <- readRDS("Output/Synonyms/01_synonyms_ias_2_or_more")

native_class <- ias_x_native_spe %>%
  filter(ias_lower %in% names_2_or_more) %>%
  distinct(binomial_iucn, Class, category, insular_endemic) %>%
  mutate(threatened = if_else(category %in% c("CR", "EN", "VU"),
                              "IAS-T","IAS-NT"))

table(native_class$Class)

table(native_class$threatened)
table(native_class$Class, native_class$threatened)

table(native_class$insular_endemic)
table(native_class$Class, native_class$insular_endemic)

table(native_class$threatened, native_class$insular_endemic)

native_name <- unique(native_class$binomial_iucn)


pbm_name <- names_2_or_more[grep("_old", names_2_or_more)]
pbm_name

names_2_or_more_clean <- gsub("_old", "", names_2_or_more)
names_2_or_more_clean <- gsub("_new", "", names_2_or_more_clean)

names_2_or_more_clean <- unique(names_2_or_more_clean)

syno_ias <- syno_ias %>% 
  distinct(accepted_name, synonym) %>%
  mutate_all(tolower)
  
# Retrieve IUCN taxonomic information for the 159 species

# iucn_search_ias <- data.frame()
# for (i in 1:length(names_2_or_more_clean)){
#   obj <- rl_search(names_2_or_more_clean[i],
#                      key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#   if(is.null(obj$result)){
#     acc_name = syno_ias$accepted_name[syno_ias$synonym==names_2_or_more_clean[i]]
#     obj_acc <- data.frame()
#     for (j in length(acc_name)){
#       obj2 <- rl_search(acc_name[j],
#                        key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")
#       obj2$result$ias_name = names_2_or_more_clean[i]
#       obj_acc <- bind_rows(obj_acc, obj2$result)
#     }
#     iucn_search_ias <- bind_rows(iucn_search_ias, obj_acc)
#   } else {
#     obj$result$ias_name = names_2_or_more_clean[i]
#     iucn_search_ias <- bind_rows(iucn_search_ias, obj$result)
#   }
#   saveRDS(iucn_search_ias,"Output/01_rl_search_ias_2_or_more")
# }

iucn_search_ias <- readRDS("Output/01_rl_search_ias_2_or_more")

# some species are not in the RedList
sum(is.na(iucn_search_ias$taxonid)) # 68
not_rl <- iucn_search_ias$ias_name[is.na(iucn_search_ias$taxonid)]

# species in the RedList
sum(table(iucn_search_ias$class)) # 94
# in which class?
table(iucn_search_ias$class)


# Are some of the 159 sp in the 100 worst alien species?

# load list of 100 worst from wikipedia
worst100 <- read.csv2("Data/100_worst_alien_species_wikipedia.csv")%>%
  mutate_all(tolower) %>%
  # correct bufo marinus
  mutate(Species = if_else(Species=="bufo marinus = rhinella marina", 
                           "rhinella marina", Species))

# for each sp of our 159, tag those in worst100
iucn_search_ias$worst100 <- "No"
iucn_search_ias$type <- "none"
for (i in 1:nrow(iucn_search_ias)){
  for (j in 1:nrow(worst100)) {
    if (iucn_search_ias$ias_name[i] == worst100$Species[j]){
      iucn_search_ias$worst100[i] <- "Yes"
      iucn_search_ias$type[i] <- worst100$Type[j]
      }
  }
}

table(iucn_search_ias$worst100) # 44 sp in worst 100
table(iucn_search_ias$type) # 13 mammals, 7 insects, & other types 

table(iucn_search_ias %>% 
        filter(is.na(taxonid)) %>%
        pull(worst100)) # 22 from the worst are not in RedList

table(iucn_search_ias %>% 
        filter(is.na(taxonid)) %>%
        pull(type)) # mostly insects





library(rgbif)

name_usage(name = names_2_or_more[1])


library(taxize)
taxize::classification(names_2_or_more[1], db = 'itis')


ritis::search_scientific("Quercus douglasii")

# chaeck for gbif data
library(rgbif)
df_count <- data.frame(ias_name = character(162), nb_occ = numeric(162))
for (i in 1:length(names_2_or_more)){
  df_count$ias_name[i] <- names_2_or_more[i]
  df_count$nb_occ[i] <- occ_search(scientificName = names_2_or_more[i])$meta$count
  print(i)
}

hist(log(df_count$nb_occ + 1), breaks = 20)

wrong_names <- df_count$ias_name[df_count$nb_occ==0]
correct_names <- c("senegalia catechu", "apis mellifera subsp. scutellata",
                   "Canis lupus subsp. dingo", "cervus elaphus", 
                   "Gallus gallus f. domesticus", "herpestes javanicus",
                   "mustela nivalis", "rhinella marina")
names(correct_names) <- wrong_names
correct_names



# sp IAS that threatened 2-3 classes
#  select only severe interactions 
ias_x_native_spe_signif <- ias_x_native_spe %>% 
  filter(severity %in% severe_inter) 

ias_class_severe <- as.data.frame(
  table(ias_x_native_spe_signif$ias_lower, ias_x_native_spe_signif$Class)) %>%
  pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>%
  rename(ias_name = Var1) %>%
  mutate(Total_inter = Aves + Mammalia + Reptilia) %>%
  mutate(Specificity = case_when(
    Aves == 0 & Mammalia == 0 ~ "Only_rept",
    Reptilia == 0 & Mammalia == 0 ~ "Only_Aves",
    Aves == 0 & Reptilia == 0 ~ "Only_mam",
    Aves > 0 & Mammalia > 0 ~ "Aves_mam",
    Reptilia > 0 & Mammalia > 0 ~ "Mam_rept",
    Aves > 0 & Reptilia > 0 ~ "Aves_rept")) %>%
  mutate(Specificity = if_else(Aves > 0 & Reptilia > 0 & Mammalia > 0,
                               "All_classes", Specificity))

table(ias_class_severe$Specificity)

table(ias_class_severe$Specificity[ias_class_severe$Total_inter == 1])

length( ias_class_severe$Total_inter[ias_class_severe$Total_inter > 1] )
sum( ias_class_severe$Total_inter[ias_class_severe$Total_inter > 1] )

names_2_or_more_signif <- unique(as.character(
  ias_class_severe$ias_name[ias_class_severe$Total_inter > 1]))
  
native_class_signif <- ias_x_native_spe_signif %>%
  filter(ias_lower %in% names_2_or_more_signif) %>%
  distinct(binomial_iucn, Class, category, insular_endemic) %>%
  mutate(threatened = if_else(category %in% c("CR", "EN", "VU"),
                              "IAS-T","IAS-NT"))

table(native_class_signif$Class)

table(native_class_signif$threatened)
table(native_class_signif$Class, native_class_signif$threatened)


table(native_class_signif$insular_endemic)
table(native_class_signif$Class, native_class_signif$insular_endemic)

table(native_class_signif$threatened, native_class_signif$insular_endemic)

native_name <- unique(native_class_signif$binomial_iucn)

# save native species list to define groups of IAS-A

# all interaction with a named IAS 
native_name_all <- unique(ias_x_native_spe %>% pull(binomial_iucn))
# Native sp associated with an IAS that interact with 2 natives or more
native_name_2more <- unique(ias_x_native_spe %>% 
                              filter(ias_lower %in% names_2_or_more) %>%
                              pull(binomial_iucn))
# Native sp that have a severe interaction with an IAS
native_name_severe <- unique(ias_x_native_spe_signif %>% pull(binomial_iucn))
# NAtive sp that have a severe interaction with an IAS that interact with 2 natives or more
native_name_severe_2more <- unique(ias_x_native_spe_signif %>%
                                     filter(ias_lower %in% names_2_or_more_signif) %>%
                                     pull(binomial_iucn))

list_native_names <- list(
  all = native_name_all,
  all2more = native_name_2more,
  severe = native_name_severe,
  severe2more = native_name_severe_2more
)

saveRDS(list_native_names, "Output/Data_clean/01_native_names_to_model")

#### Test open Hof file

# install.packages('ncdf4')
library(raster)
library(ncdf4)
r1 <- raster("Data/bioscen1.5-sdm-gam_gfdl-esm2m_ewembi_historical_nosoc_co2_thrmammalsumprob_global_30year-mean_1990_1990.nc")
r2 <- raster("Data/bioscen1.5-sdm-gam_gfdl-esm2m_ewembi_rcp26_nosoc_co2_thrmammalsumprob_global_30year-mean_2009_2080.nc")
r3 
plot(r2)
str(r2)
r2@data

i = "Data/bioscen1.5-sdm-gam_gfdl-esm2m_ewembi_historical_nosoc_co2_thrmammalsumprob_global_30year-mean_1990_1990.nc"
nc <- nc_open(i)
v <- nc$var[[nc$nvars]]
test <- stack(i, varname=v$name)
test2 = stack(i, varname = "Acanthixalus_spinosus")
print(nlayers(test))
plot(test[[sample(1:nlayers(test),1)]])
