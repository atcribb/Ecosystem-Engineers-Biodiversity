#Stratigraphic binning to make columns for stratigraphic binng 
#             column names = $stage
#.                           #series
#Author: Alison Cribb (A.T.Cribb@soton.ac.uk)
#Created: 20 September 2023
#Last edited: 1 November 2023

library(divDyn)
library(dplyr)

data(stages)
#colnames(stages)[colnames(stages)=="stage"] <- "name"
data(keys)

stage_binning <- function(pbdb_data){
  
  #Following ddPhanero (Kocscis et al.)
  stgMin <- categorize(pbdb_data[,'early_interval'], keys$stgInt)
  stgMin <- as.numeric(stgMin)
  stgMax <- categorize(pbdb_data[,'late_interval'], keys$stgInt)
  stgMax <- as.numeric(stgMax)
  
  stgCondition <- c(
    which(stgMax==stgMin),
    which(stgMax==-1))
  pbdb_data$stg[stgCondition] <- stgMin[stgCondition]
  
  return(pbdb_data)
  
}



