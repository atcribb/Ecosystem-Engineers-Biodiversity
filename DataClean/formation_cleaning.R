#Author: Alison Cribb
#Date created: 6 March 2024
#Last edited:
#Summary: Cleans PBDB data at the formation level to remove synonyms, uncertain formations, etc. 

library(readr)

setwd("~/Desktop/Manucripts/EcosystemEngineering_Biodiversity/Data")

load('Phanerozoic_marine_cleaned_binned.RData')
all_data <- marine_cleaned_binned
formation_cleaning <- read.csv('formation_sorting.csv')

colnames(all_data)[which(colnames(all_data)=='formation')] <- 'old_formation'
all_data <- subset(all_data, !(old_formation==''))

all_data$formation <- rep(NA, nrow(all_data))
for(i in 1:nrow(all_data)){
  this.old.formation <- all_data$old_formation[i]
  corrected.formation.row <- subset(formation_cleaning, old_formation==this.old.formation)
  corrected.formation <- corrected.formation.row$change_to
  if(length(corrected.formation>0)){
    all_data$formation[i] <- corrected.formation[1]
  }
}
save(all_data, file='Phanerozoic_clean_final.RData')
