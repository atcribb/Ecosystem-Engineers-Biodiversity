#Creates reference dataset for bioturbating ecosystem engineers

#=== load data ===#
load('Data/Phanerozoic_clean_final.RData')

#===== Find fossil assemblages with and without bioturbators ====#
table(all_data$life_habit) #use life_habit to find infauna 
infaunal_lifehabits <- c('deep infaunal', 'infaunal', 'semi-infaunal', 'semi-infaunal, solitary', 'shallow-infaunal',
                         'infaunal, depth=deep', 'semi-infaunal, gregarious', 'semi-infaunal, solitary',
                         'shallow infaunal')
infauna_data <- subset(all_data, life_habit %in% infaunal_lifehabits)
infauna_EE_genera <- unique(infauna_data$genus)

#=== and also add epifaunal sediment bulldozers  ===#
epifaunal_data <- subset(all_data, life_habit=='epifaunal')
epifaunal_data_mobile <- subset(epifaunal_data, motility=='actively mobile') #conservatively, just the actively mobile epifauna
EE_epifauna_diets <- c('detritivore, grazer', 'grazer', 'grazer, deposit feeder') #conservative guesses for what's eating the top surface of the seafloor
bulldozers_EE_data <- subset(epifaunal_data_mobile, diet %in% EE_epifauna_diets)
bulldozers_EE_genera <- unique(epifauna_EE_data$genus)

#=== combine for final list of benthic->infauna ecosystem engineer genera and formations ===# 
ecoeng_genera <- unique(c(infauna_EE_genera, bulldozers_EE_genera))
initial_grab_data <- subset(all_data, genus %in% ecoeng_genera)
table(initial_grab_data$phylum) #some of these we can tell shouldn't be included
keep_phyla <- c('Annelida', 'Arthropoda', 'Brachiopoda', 'Echinodermata', 'Mollusca')

#final final data to use 
bioturbators_data <- subset(initial_grab_data, phylum %in% keep_phyla)
length(unique(bioturbators_data$genus))
length(unique(
  subset(bioturbators_data, !is.na(formation))$formation
))

save(bioturbators_data, file='Data/Bioturbators_data.RData')