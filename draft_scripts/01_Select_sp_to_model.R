# Focus on species associated to native amniotes
# Count IAS associated to each group/all groups

rm(list=ls())

library(tidyr)
library(dplyr)
library(rredlist)
library(tidyverse)
library(giscoR)
library(rgbif)
library(taxize)


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
saveRDS(names_2_or_more, "Output/01_names_2_or_more_162")
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
saveRDS(names_2_or_more_clean, "Output/01_names_2_or_more_clean")


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


#-------------- Are some of the 159 sp in the 100 worst alien species?

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


#-------------- Are some of the 159 in Nguyen & Leung 2022 ?

sdm647 <- read.csv2("Data/sdm_ias_647_N&L2020.csv")
sdm647$species_lower <- tolower(sdm647$Species)

sp_in_sdm <- names_2_or_more_clean[names_2_or_more_clean %in% sdm647$species_lower]
sp_not_in_sdm <- names_2_or_more_clean[!(names_2_or_more_clean %in% sdm647$species_lower)]


# for each sp of our 159, tag those in sdm647
iucn_search_ias$sdm647 <- "No"
iucn_search_ias$taxa <- "na"
for (i in 1:nrow(iucn_search_ias)){
  for (j in 1:nrow(sdm647)) {
    if (iucn_search_ias$ias_name[i] == sdm647$species_lower[j]){
      iucn_search_ias$sdm647[i] <- "Yes"
      iucn_search_ias$taxa[i] <- sdm647$Taxa[j]
    }
  }
}

table(iucn_search_ias$sdm647)
table(iucn_search_ias$sdm647, iucn_search_ias$worst100)

# ---------------- Find taxonomy of all 159 species

# taxo_name_ias <- data.frame()
# for (i in names_2_or_more_clean){
#   taxo_name_ias <- bind_rows(
#     taxo_name_ias,
#     taxize::tax_name(i, get = c('kingdom','class','order',"family"), 
#                      db='ncbi'))
# }
# 
# saveRDS(taxo_name_ias,"Output/01_name_ias_2_or_more_taxo_ncbi")

taxo_name_ias <- readRDS("Output/01_name_ias_2_or_more_taxo_ncbi")

#taxize::rank_ref

hist(table(taxo_name_ias$class))

agg_class <- taxo_name_ias %>% 
  group_by(class) %>%
  summarise(count=n())

ggplot(agg_class, aes(x = reorder(class, count), y = count)) +
  geom_bar(stat = 'identity') + xlab("")+
  theme_classic() + coord_flip()

# agg_order <- taxo_name_ias %>% 
#   group_by(order) %>%
#   summarise(count=n())
# ggplot(agg_order, aes(x = reorder(order, count), y = count)) +
#   geom_bar(stat = 'identity') + xlab("")+
#   theme_classic() + coord_flip()
# 
# 

taxo_name_ias_in_db <- left_join(
  taxo_name_ias %>% rename(binomial = query),
  iucn_search_ias %>% select(ias_name, worst100, sdm647) %>%
    rename(binomial = ias_name),
  by = "binomial")

write.table(taxo_name_ias, "Output/01_159_IAS_interact_2_or_more_tetrap.txt")


# chaeck for gbif data
library(rgbif)
# df_count <- taxo_name_ias_in_db %>% 
#   mutate(nb_occ = numeric(nrow(taxo_name_ias_in_db)))
# for (i in 1:nrow(df_count)){
#   df_count$nb_occ[i] <- occ_search(scientificName = df_count$binomial[i])$meta$count
#   print(i)
# }
# 
# hist(log(df_count$nb_occ + 1), breaks = 20)
# 
# wrong_names <- df_count$binomial[df_count$nb_occ==0]
# 
# correct_names <- c("acacia catechu"="senegalia catechu", 
#                    "apis mellifera ssp. scutellata" = "Apis mellifera subsp. scutellata",
#                    "canis lupus ssp. dingo" ="Canis lupus subsp. dingo", 
#                    "gallus gallus ssp. domesticus" = "Gallus gallus f. domesticus")
# correct_names
# 
# for (i in 1:nrow(df_count)){
#   for (j in 1:length(correct_names)){
#     if (df_count$binomial[i] == names(correct_names)[j]){
#       df_count$nb_occ[i] <- occ_search(scientificName = correct_names[j])$meta$count
#       print(correct_names[j])
#     }
#   }
# }
# 
# saveRDS(df_count, "Output/01_159_IAS_count_gbif")

df_count <- readRDS("Output/01_159_IAS_count_gbif")

hist(log(df_count$nb_occ + 1), breaks = 20)


n_occ_100 <- read.csv2("Data/100_worst_CB_list_2013_N_occ.csv")


head(df_count)

df_count_rl <- df_count %>% 
  mutate(redlist = if_else(
    binomial %in% pull(iucn_search_ias %>% filter(!is.na(taxonid)), ias_name), 
    "Yes", "No"))


table(df_count_rl$redlist)
table(df_count_rl$worst100)
table(pull(df_count_rl %>% filter(worst100=="No"), sdm647))
table(pull(df_count_rl %>% filter(worst100=="No"), redlist))
table(pull(df_count_rl %>% filter(worst100=="No" & sdm647=="No"), redlist))
# 30 species are not either in the redlist, in worst 100 nor sdm647

pull(df_count_rl %>% filter(worst100=="No" & sdm647=="No" & redlist == "No"), binomial)

# look if those species are in gisd or cabi:
cabi <- c("andropogon gayanus", "apis mellifera", "apis mellifera ssp. scutellata",
          "bacillus anthracis", "bos taurus", "bubalus bubalis", "camelus dromedarius",
          "canis familiaris", "ctenopharyngodon idella", "cyathea cooperi",
          "dolichandra unguis-cati", "equus asinus", "equus caballus", 
          "mustela furo", "ovis aries", "pasteurella multocida", "philornis downsi",
          "rubus argutus", "schistocerca nitens", "solenopsis geminata", 
          "sophonia rufofascia", "themeda quadrivalvis", "toxoplasma gondii",
          "yersinia pestis")

gisd <- c("apis mellifera ssp. scutellata", "bos taurus", 
          "bubalus bubalis", "camelus dromedarius", "canis familiaris",
          "ctenopharyngodon idella", "cyathea cooperi", "equus asinus",
          "equus caballus", "gallus gallus ssp. domesticus", "mustela furo",
          "ovis aries", "pasteurella multocida", "philornis downsi", 
          "solenopsis geminata", "yersinia pestis")

# COMMENTS
# canis familiaris => native area = australia??
# canis lupus ssp. dingo = canis lupus in gisd
# gallus gallus ssp. domesticus = gallus gallus in gisd
# sophonia rufofascia is sophonia orientalis in cabi
# toxoplasma gondii & yersinia pestis have a cabi/gisd sheet but not with native range

other_ref <- c("brachiaria decumbens", "canis lupus ssp. dingo",
               "sophonia rufofascia")
# https://www.feedipedia.org/node/489
# https://animaldiversity.org/site/accounts/information/Canis_lupus_dingo.html
# https://www.researchgate.net/publication/233680334_
# Effect_of_Sophonia_rufofascia_Homoptera_Cicadellidae_on_Guava_Production_in_Hawaii

none <- c("Clostridium botulinum", "erysipelothrix rhusiopathiae", "sus domesticus",
          "toxoplasma gondii")
# Clostridium botulinum = bacteria, see on gbif if "native" area?
# erysipelothrix rhusiopathiae = bacteria
# sus domesticus to bind with sus scrofa?

data_72 <- read.delim("Data/data72R.txt")
colnames(data_72)
str(data_72)

# select only global studies
table(data_72$scale2)
table(data_72$scale3)

table(pull(data_72 %>%
             filter(scale3 == "world"),
           taxonomic_group))

sp <- unique(pull(data_72 %>%
             filter(scale3 == "world"),
           invasive_species))

df_count_all <- df_count_rl %>%
  mutate(in_72_articles = if_else(binomial %in% tolower(sp), "Yes", "No"),
         in_CABI = if_else(binomial %in% cabi, "Yes", "No"),
         in_gisd = if_else(binomial %in% gisd, "Yes", "No"),
         in_other = if_else(binomial %in% other_ref, "Yes", "No")) %>%
  mutate(nada = if_else(redlist =="No" & worst100 == "No" & sdm647 == "No" & 
                          in_72_articles == "No" &in_CABI == "No" & 
                          in_gisd == "No" & in_other == "No",
                        "nada", "ok"))
table(df_count_all$nada)

colnames(df_count_all)

df_count_all %>% filter(nada=="nada") %>% pull(binomial)

# vu avec Céline
# keep only species that are terrestrial (remove freshwater fishes)
# remove bacteria and diseases - keep parasites?
# remove sus domesticus but change it to sus scrofa in ias_x_native

df_count_terr <- df_count_rl %>%
  # add a class for 
  mutate(kingdom = if_else(
    binomial %in% c("apis mellifera ssp. scutellata", "canis lupus ssp. dingo",
                    "gallus gallus ssp. domesticus"), 
    "Metazoa", kingdom),
    class = if_else(binomial == "apis mellifera ssp. scutellata", "Insecta", class),
    class = if_else(binomial == "canis lupus ssp. dingo", "Mammalia", class),
    class = if_else(binomial == "gallus gallus ssp. domesticus","Aves", class)) %>%
  # remove species that have no info nowhere
  # or that are bacteria / protista / parasite
  filter(!(is.na(kingdom))) %>%
  # remove fishes
  filter(class != "Actinopteri") %>%
  filter(binomial != "sus domesticus")


setdiff( df_count_rl$binomial, df_count_terr$binomial)

agg_class <- df_count_terr %>%
  group_by(class) %>% 
  summarise(count= n())
ggplot(agg_class, aes(x = reorder(class, count), y = count)) +
  geom_bar(stat = 'identity') + xlab("")+
  theme_classic() + coord_flip()


# get their gbif ID
sp_list_for_gbif_ok <- df_count_terr$binomial
correct_names <- c("acacia catechu"="senegalia catechu", 
                   "apis mellifera ssp. scutellata" = "Apis mellifera subsp. scutellata",
                   "canis lupus ssp. dingo" ="Canis lupus subsp. dingo",
                   "gallus gallus ssp. domesticus" = "Gallus gallus f. domesticus")
sp_list_for_gbif_ok[sp_list_for_gbif_ok %in% names(correct_names)] <- 
  correct_names

df_count_terr$binomial_corr <- sp_list_for_gbif_ok
# test <- c("senegalia catechu", "sus scrofa", "rubus niveus")


gbif_taxon_keys_terr <-
  df_count_terr %>%
  #  filter(binomial_corr %in% test) %>% # use fewer names if you want to just test
  pull(binomial_corr) %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %>% # combine all data.frames into one
  filter(status != "DOUBTFUL" & matchtype != "FUZZY")
# choose to download all data from gbif that coul match our species names
# then make a cleaning to keep only occurrences linked to each sp
# %>% filter(matchtype == "EXACT" & status == "ACCEPTED") # get only accepted and matched names

length(unique(gbif_taxon_keys_terr$original_sciname))
length(unique(gbif_taxon_keys_terr$canonicalname))

nrow(gbif_taxon_keys_terr %>% distinct(original_sciname, order))


gbif_taxon_keys_terr %>%
  distinct(canonicalname, original_sciname)

gbif_taxon_keys_terr %>% 
  distinct(original_sciname, order)

gbif_taxon_keys_terr_hb <- gbif_taxon_keys_terr %>%
  filter(!(species %in% c("Vachellia nilotica", "Podalonia hirsuta",
                          "Cenchrus brownii", "Miconia hirta", "Dalbergia rostrata",
                          "Mikania cynanchifolia", "Rubus pedunculosus",
                          "Salix pentandra"
                          ))) %>%
  filter(!(scientificname %in% c("Bothriochloa pertusa (L.) Maire", 
                                 "Cenchrus echinatus Cav.", 
                                 "Gallus gallus domesticus (Linnaeus, 1758)",
                                 "Hedychium gardnerianum Wall.",
                                 "Hedychium gardnerianum Sheppard",
                                 "Hemidactylus frenatus Boie In Schlegel",
                                 "Hemidactylus frenatus Kuhl & Van Hasselt In Schlegel",
                                 "Leucaena leucocephala (Lam.) Dewit",
                                 "Lycodon capucinus Boie In F.Boie, 1827",
                                 "Myrica faya Dryand.",
                                 "Psidium cattleianum Afzel.",
                                 "Viverricula indica (Desmarest, 1804)")))

length(unique(gbif_taxon_keys_terr_hb$original_sciname))
length(unique(gbif_taxon_keys_terr_hb$canonicalname))

nrow(gbif_taxon_keys_terr_hb %>% distinct(original_sciname, order))

write.csv2(gbif_taxon_keys_terr_hb, "Output/Data_clean/01_145_IAS_to_model.csv",
           row.names = F)


##################### SP IN DASCO ?? #######################


ias145 <- read.csv2("Output/Data_clean/01_145_IAS_to_model.csv")

sp_dasco <- read.csv2("Z:/THESE/5_Data/Alien_data/Dasco/DASCO_AlienRegions_SInAS_2.4.1.csv")

dsc <- sp_dasco %>% distinct(taxon, scientificName)
unique(sp_dasco$location)

in_dasco <- unique(sp_dasco %>% 
  filter(taxon %in% ias145$canonicalname) %>%
  pull(taxon))

# what are the species absent from dasco?
setdiff(ias145$canonicalname, in_dasco)

# canis lupus familiaris = canis familiaris
sp_dasco$taxon[grep("canis", tolower(sp_dasco$taxon))]
# herpestes javanicus auropunctatus = herpestes auropunctatus
sp_dasco$taxon[grep("herpestes", tolower(sp_dasco$taxon))]


#try to find species name as in gbif ?
sp_no_dasco <- ias145 %>%
   filter(canonicalname %in% setdiff(ias145$canonicalname, in_dasco))

sp_dasco %>% filter(scientificName %in% sp_no_dasco$scientificname)


in_dasco <- c(in_dasco, "Canis lupus familiaris", "Herpestes javanicus auropunctatus")


occ_dasc_obis <- read.csv2("Z:/THESE/5_Data/Alien_data/Dasco/DASCO_OBISCoords_SInAS_2.4.1.gz", sep = ",")
occ_dasc_gbif <- read.csv2("Z:/THESE/5_Data/Alien_data/Dasco/DASCO_GBIFCoords_SInAS_2.4.1.gz", sep = ",")



occ_dasc_obis_sp <- occ_dasc_obis %>%
  filter(taxon %in% ias145$species)

spkey <- ias145$specieskey
usagekey <- ias145$usagekey # considère moins d'espèces 

occ_dasc_gbif_sp <- occ_dasc_gbif %>%
  filter(speciesKey %in% spkey)

length(unique(occ_dasc_gbif_sp$speciesKey))

unique(occ_dasc_obis_sp$taxon)


# bind obis and gbif records
sp_key_taxon <- ias145 %>%
  distinct(specieskey, species) %>%
  rename(speciesKey = specieskey,
         taxon = species)
# on tombe ici à 140 espèces uniques, parce que parmi nos 145 il y a des sous espèces


colnames(occ_dasc_gbif_sp)
colnames(occ_dasc_obis_sp)

occ_dasc_gbif_sp_col <- left_join(
  occ_dasc_gbif_sp, sp_key_taxon,
  by = "speciesKey"
)
occ_dasc_obis_sp_col <- left_join(
  occ_dasc_obis_sp, sp_key_taxon,
  by = "taxon"
)

occ_all <- bind_rows(occ_dasc_gbif_sp_col, occ_dasc_obis_sp_col) %>%
  mutate(LONG = as.numeric(Longitude),
         LAT = as.numeric(Latitude)) %>%
  select(-c(Latitude, Longitude))
length(unique(occ_all$taxon))

nb_occ_inv_tax <- table(occ_all$taxon)
hist(log(table(occ_all$taxon)))

library(sp)
library(raster)

convertcoord <- function(x){
  WGScoor <-  x
  sp::coordinates(WGScoor)=~LONG+LAT
  sp::proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
  return(WGScoor)}

ias_spatial_129 <- convertcoord(occ_all)

saveRDS(ias_spatial_129, "Output/focus_ias_range/01_ias129_occ_dasco_alien_range")


# ----------- From points to raster for each sp 

ias_spatial_129 <- readRDS("Output/focus_ias_range/01_ias129_occ_dasco_alien_range")

# world map to filter out oceanic ranges ?
# library(rnaturalearth)
# worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sp")

# INPUT : a spatial point object and a raster filled with 0
# OUTPUT : a df containing all the cells with a point
extract_cells_pts <- function(pts, rast){
  
  # extract cells with a point inside
  c_pts <- extract(rast, pts, cellnumbers=TRUE)[,"cells"]
  
  # convert raster cells with values to 1 
  rast_val <- rast
  rast_val[c_pts] <- 1

  # for saving memory, save a df containing all cells with values
  df_val <- as.data.frame(rast_val, xy = TRUE) %>%
    # remove empty cells
    dplyr::filter(layer == 1) %>%
    dplyr::mutate(binomial = unique(pts$taxon))
  
  return(df_val)
}
# herp = subset(ias_spatial_129, taxon == "Herpestes javanicus")
# test <- extract_cells_pts(herp, raster0.1)

raster0.1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90,
                    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    resolution = 0.1, vals=NULL)
raster0.1[] <- 0

library(pbapply)
taxon <- as.list(unique(ias_spatial_129$taxon))

df_129_ias <- pblapply(taxon, function(x){
  pts_sp <- subset(ias_spatial_129, taxon == x)
  cells <- extract_cells_pts(pts_sp, raster0.1)
  return(cells)
}) # takes 02m 13s

saveRDS(df_129_ias, "Output/focus_ias_range/01_ias129_raster_df")

#----------- MAP IAS SPECIES -----------

df_129_ias <- readRDS("Output/focus_ias_range/01_ias129_raster_df")

df_all <- bind_rows(df_129_ias) 


# SR of 129 alien species
ias_all_agg <- df_all %>%
  group_by(x, y) %>%
  summarise(SR_tot = n())


SR_ias<- ggplot(data = ias_all_agg) +
  geom_raster(aes(x = x, y = y, fill = SR_tot)) +
  scale_fill_gradient(
    low = "yellow", 
    high = "red")+
  #geom_sf(data = worldMap, alpha = 0.1) +
  theme_classic() +
  labs(title = "Species richness of the target IAS",
       subtitle = paste0("Total number of IAS with at least one pixel: ",
                         length(unique(df_all$binomial))),
       x = "Longitude", y = "Latitude")
SR_ias




##### from df to raster to extract SR values for polygon islands
library(raster)
spg <- ias_all_agg
coordinates(spg) <- ~ x + y
# coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE
# coerce to raster
ras_ias <- raster(spg)

# world map to filter out oceanic ranges
library(rnaturalearth)
worldMap <- ne_countries(scale = "medium", type = "countries", returnclass = "sp")

# pour chaque entité, extract the cells that fall into
ext_world <- extract(ras_ias, worldMap)
names(ext_world) <- worldMap$name

# takes max
ext_world_max <- lapply(ext_world, function(x){max(na.omit(x))})

ext_world_max_df <- as.data.frame(unlist(ext_world_max)) %>%
  mutate(id = worldMap$name)

names(ext_world_max_df) <- c("max_ias", "id")

worldMap.fort <- fortify(worldMap, region='name') 
worldMap_df <- left_join(worldMap.fort, ext_world_max_df, by = "id")


ggplot(worldMap_df)+
  geom_polygon(aes(x=long, y=lat, group=group, fill = max_ias)) +
  scale_fill_gradient(low = "white", high = "#AE123A")




sp_no_dasco <- setdiff(sp_key_taxon$taxon, unique(df_all$binomial))
sp_no_dasco



occ_countries <- data.frame()
for (name in sp_no_dasco){
    obj <- rl_occ_country(name,
                          key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")$result
    obj$binomial = name
    occ_countries <- bind_rows(occ_countries, obj)
}
length(unique(occ_countries$binomial))

saveRDS(occ_countries,"Output/01_countries_occ_ias_iucn")


# For species in IUCN, compute actual exposure
# all grid cells were an invasive sp occur outside of its native range (IUCN)
# or within its invasive range (GISD/ CABI)

# occ_countries <- data.frame()
# for (i in 1:nrow(iucn_search_ias)){
#   name <- iucn_search_ias$scientific_name[i]
#   if (!is.na(name)){
#     obj <- rl_occ_country(name,
#                           key = "0e9cc2da03be72f04b5ddb679f769bc29834a110579ccdbeb54b79c22d3dd53d")$result
#     obj$binomial = name
#     occ_countries <- bind_rows(occ_countries, obj)
#   }
# }
# saveRDS(occ_countries,"Output/01_countries_occ_ias_iucn")

occ_countries <- readRDS("Output/01_countries_occ_ias_iucn")
native_countries <- occ_countries %>%
  filter(origin %in% c("Native","Reintroduced","Vagrant"))

# for each sp, collect all countries from native area
native_area_sp <- vector(mode = "list", 
                         length = length(unique(native_countries$binomial)) )
names(native_area_sp) <- unique(native_countries$binomial)
# for each sp, collect all occurrences in gbif
gbif_occur <- native_area_sp


# for(i in unique(native_countries$binomial)){
#   native_area_sp[[i]] <- gisco_get_countries(
#     country = native_countries %>%
#       filter(binomial==i) %>% pull(code))
# 
#     gbif_occur[[i]] <- occ_search(scientificName = i)
#     
#     print(i)
# }
# native_area_sp$`Senegalia catechu`


# Retrieve and save GBIF data for all 159 sp to get actual distribution
sp_list_for_gbif <- df_count$binomial
correct_names <- c("acacia catechu"="senegalia catechu", 
                   "apis mellifera ssp. scutellata" = "Apis mellifera subsp. scutellata",
                   "canis lupus ssp. dingo" ="Canis lupus subsp. dingo",
                   "gallus gallus ssp. domesticus" = "Gallus gallus f. domesticus")
sp_list_for_gbif[sp_list_for_gbif %in% names(correct_names)] <- 
  correct_names

df_count_rl$binomial_corr <- sp_list_for_gbif
# test <- c("senegalia catechu", "sus scrofa", "rubus niveus")

# fill in your gbif.org credentials 
user <- "claramarino" # your gbif.org username 
pwd <- "data2021" # your gbif.org password
email <- "claramarino665@gmail.com" # your email 

library(taxize)


gbif_taxon_keys <-
  df_count_rl %>%
#  filter(binomial_corr %in% test) %>% # use fewer names if you want to just test
  pull(binomial_corr) %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() # combine all data.frames into one
# choose to download all data from gbif that coul match our species names
# then make a cleaning to keep only occurrences linked to each sp
# %>% filter(matchtype == "EXACT" & status == "ACCEPTED") # get only accepted and matched names

length(unique(gbif_taxon_keys$original_sciname))
length(unique(gbif_taxon_keys$canonicalname))

nrow(gbif_taxon_keys %>% distinct(original_sciname, order))

agg_class <- gbif_taxon_keys %>% distinct(original_sciname, phylum) %>% 
  group_by(phylum) %>%
  summarise(count=n())

ggplot(agg_class, aes(x = reorder(phylum, count), y = count)) +
  geom_bar(stat = 'identity') + xlab("")+
  theme_classic() + coord_flip()



occ_download(
  pred_in("taxonKey", unique(gbif_taxon_keys$usagekey)),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)



# Open very large file with R
occ_3_sp <- readr::read_tsv("Data/GBIF_occurrences/0228384-210914110416597.csv")
occ_159_sp <- readr::read_tsv("Data/GBIF_occurrences/0228550-210914110416597.csv")

# cannot be open with r because too heavy => to handle with readr::read_tsv_chunked()




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
