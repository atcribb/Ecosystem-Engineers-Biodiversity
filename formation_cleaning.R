#Author: Alison Cribb
#Date created: 6 March 2024
#Last edited:
#Summary: Cleans PBDB data at the formation level to remove synonyms, uncertain formations, etc. 

library(readr)


setwd("~/Desktop/Manucripts/EcosystemEngineering_Biodiversity/Data")

load('Phanerozoic_marine_cleaned_binned.RData')
load('Reef_Ecosystem_Engineers.RData')
all_data <- marine_cleaned_binned
formation_cleaning <- read.csv('formation_sorting.csv')

#as an example of why we have to do this...
homerian <- subset(all_data, stage=='Homerian')
table(homerian$formation)

colnames(all_reef_builders)[which(colnames(all_reef_builders)=='formation')] <- 'old_formation'
all_reef_builders <- subset(all_reef_builders, !(old_formation==''))

all_reef_builders$formation <- rep(NA, nrow(all_reef_builders))

for(i in 1:nrow(all_reef_builders)){
  
  this.old.formation <- all_reef_builders$old_formation[i]
  corrected.formation <- subset(formation_cleaning, old_formation==this.old.formation)$change_to
  if(length(corrected.formation>0)){
  all_reef_builders$formation[i] <- corrected.formation[1]
  }
}
nrow(subset(all_reef_builders, is.na(formation)))
save(all_reef_builders, file='Reef_Ecosystem_Engineers_final.RData')


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
