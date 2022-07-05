# RISK find alien range from GRIIS using Anna's script

######################## IASNH interactions ############################
########################################################################
## script by Anna J. Turbelin
## created on 05 October 2021

##package for GBIF data cleaning
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_
library(gdistance)
library(spThin)
library(rworldmap)
library(rworldxtra)
library(maptools)
library(maps)
library(data.table)
library(bRacatus)
library(shinyjs)
library(stringr)
library(countrycode)
library(CoordinateCleaner)
library(readxl)
library(sp)
library(rgdal)
library(raster)
library(dismo)
library(rJava)

options(stringsAsFactors=FALSE)


rm(list=ls())

# load species list
gbif_taxo_2 <- readRDS("Output/Taxonomy/RISK_01_taxo_gbif_2")


SpList <- unique(gbif_taxo_2$species) ##list of alien species to search for


## check alien range
path_griis <- "Z:/THESE/5_Data/Alien_data/GRIIS_data_Anna"
list_dir <- list.dirs(path_griis, recursive=FALSE) #list folders' names in GRIIS_data folder - contains all GRIIS datasets except protected areas  


df <- data.frame(matrix(nrow = 0, ncol= 6 ), stringsAsFactors = FALSE) ## create blank dataframe
colnames(df) <- c("scientificName","locationID","countryCode", "occurrenceStatus", "establishmentMeans","Species_or") ## name column of blank df "Species_or" is the original name from the species list  


for (j in 1:length(SpList)){
  
  sp <- SpList[j] 
  
  for (i in 1:length(list_dir)){
    
    if (j == 1){l <- i} else {l <- i+ ((j-1)*221)} #inventory list number so that the rows do not overlap
    
    GRIIS_file_loc <- list_dir[i]
    GRIIS_file <- list.files(paste0(GRIIS_file_loc), recursive=FALSE)
    GRIIS_file <- GRIIS_file[grepl(".csv", GRIIS_file)] #get list of files that end with .csv
    GRIIS_df <- fread(paste0(GRIIS_file_loc,"/",GRIIS_file ),sep = ",") #read file
    
    sp_row <- which(grepl(sp,GRIIS_df$scientificName))
    sp_row1 <- which(grepl(sp,GRIIS_df$acceptedNameUsage))
    scientificName <- which(colnames(GRIIS_df)=="scientificName")
    locationID <- which(colnames(GRIIS_df)=="locationID")
    countryCode <-which(colnames(GRIIS_df)=="countryCode")
    occurrenceStatus <-which(colnames(GRIIS_df)=="occurrenceStatus")
    establishmentMeans <-which(colnames(GRIIS_df)=="establishmentMeans")
    
    if(length(sp_row)>0){
      if(length(sp_row)==1){
        df[l,1] <- GRIIS_df[sp_row,scientificName]
        if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row,locationID]}
        df[l,3] <- GRIIS_df[sp_row,countryCode]
        if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row,occurrenceStatus]}else{df[l,4] <- "NA"}
        df[l,5] <- GRIIS_df[sp_row,establishmentMeans]
        df[l,6] <- sp
      }else
        if(length(sp_row)>1){
          # sp_row2 <- which(GRIIS_df$taxonRank == "SPECIES" & grepl(sp,GRIIS_df$scientificName))
          sp_row2 <- min(sp_row)
          # if(length(sp_row2==0)){sp_row2 <- min(sp_row)}
          df[l,1] <- GRIIS_df[sp_row2,scientificName]
          if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row2,locationID]}
          df[l,3] <- GRIIS_df[sp_row2,countryCode]
          if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row2,occurrenceStatus]}else{df[l,4] <- "NA"}
          df[l,5] <- GRIIS_df[sp_row2,establishmentMeans]
          df[l,6] <- sp
        }
      
    }
    
    if(length(sp_row1)>0){
      if(length(sp_row1)==1){
        df[l,1] <- GRIIS_df[sp_row1,scientificName]
        if(length(locationID)>0){ df[l,2] <- GRIIS_df[sp_row1,locationID]}
        df[l,3] <- GRIIS_df[sp_row1,countryCode]
        if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row1,occurrenceStatus]}else{df[l,4] <- "NA"}
        df[l,5] <- GRIIS_df[sp_row1,establishmentMeans]
        df[l,6] <- sp
      }else
        if(length(sp_row1)>1){
          sp_row3 <- sp_row1[1]
          df[l,1] <- GRIIS_df[sp_row3,scientificName]
          if(length(locationID)>0){df[l,2] <- GRIIS_df[sp_row3,locationID]}
          df[l,3] <- GRIIS_df[sp_row3,countryCode]
          if(length(occurrenceStatus)>0){ df[l,4] <- GRIIS_df[sp_row3,occurrenceStatus]}else{df[l,4] <- "NA"}
          df[l,5] <- GRIIS_df[sp_row3,establishmentMeans]
          df[l,6] <- sp
        }
    }
    else {next}
    
  }
  
  rm(list=setdiff(ls(), c("df", "SpList", "path_griis"))) ## remove stuff from R memory
  list_dir <- list.dirs(path_griis, recursive=FALSE)
}

df_clean <- na.omit(df)

str(df_clean)

length(unique(df_clean$Species_or))

write.csv(df_clean, "IASNH_sp_invasiveRange.csv")


