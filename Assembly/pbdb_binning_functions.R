#Stratigraphic binning to make columns for stratigraphic binng 
#             column names = $stage
#.                           #series
#Author: Alison Cribb (A.T.Cribb@soton.ac.uk)
#Created: 20 September 2023
#Last edited: 1 November 2023

setwd("~/Desktop/Manucripts/Phanerozoic_climate_tracking/Functions")

library(divDyn)
library(dplyr)

data(stages)
data(keys)

stage_binning <- function(Paleozoic_data){
  
  #Following ddPhanero (Kocscis et al.)
  #for the Paleozoic_data
  stgMin <- categorize(Paleozoic_data[,'early_interval'], keys$stgInt)
  stgMin <- as.numeric(stgMin)
  stgMax <- categorize(Paleozoic_data[,'late_interval'], keys$stgInt)
  stgMax <- as.numeric(stgMax)
  
  Paleozoic_data$stg <- rep(NA, nrow(Paleozoic_data))
  
  stgCondition <- c(
    which(stgMax==stgMin),
    which(stgMax==-1))
  Paleozoic_data$stg[stgCondition] <- stgMin[stgCondition]
  
  #remove any data without a stage assignment 
  Paleozoic_data <- Paleozoic_data[!is.na(Paleozoic_data$stg),] #losing 32,282 entries... lol but sitll keeping >100k
  
  #and finally put names to stgs 
  for(i in 1:length(Paleozoic_data$stg)){
    #find the stage number assinged from earl and late interval
    stage_number <- Paleozoic_data$stg[i]
    #what is the corrseponding stage in stages?
    stagename <- stages$stage[stage_number]
    #put that back in the dataset
    Paleozoic_data$stage[i] <- stagename
  }
  
  return(Paleozoic_data)
  
}




