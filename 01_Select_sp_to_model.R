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
