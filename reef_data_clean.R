#Author: AT Cribb
#Created: 27 November 2023
#Last edited:

#Summary: create a reef builder dataset from the cleaned and binned marine phanerozoic data

#===================================================================#
setwd("~/Desktop/Manucripts/Palass_ecosystemengineering")

#=== load data ===#
load('Data/Phanerozoic_marine_cleaned_binned.RData')
all_data <- marine_cleaned_binned

#=== reef grabs ===#
#archaeocyathids 
archaeocyath_family <- c('Archaeocyathidae') #archaeocyathids 
archaeocyath_finds <- subset(all_data, all_data$family == archaeocyath_family)
archaeocyath_finds$REE_classification <- 'Archaeocyaths'

#stromatoporoids
stromatoporoid_class <- c('Stromatoporoidea') #stromatoporoids 
stromatoporoid_finds <- subset(all_data, all_data$class == stromatoporoid_class)
stromatoporoid_finds$REE_classification <- 'Stromatoporoids'

#glass sponge reefs
hexactinellida_class <- c('Hexactinellida') #glass sponge reefs
hexactinellida_finds <- subset(all_data, all_data$class == hexactinellida_class)
hexactinellida_finds$REE_classification <- 'Glass_sponges'

#rudist bivalves
rudist_fams <- c("Caprinulidae","Caprotinidae","Diceratidae","Hippuritidae","Monopleuridae","Plagioptychidae","Polyconitidae","Radiolitidae","Trechmannellidae") #rudist corals 
rudist_finds <- subset(all_data, all_data$family %in% rudist_fams)
rudist_finds$REE_classification <- 'Rudist_bivalves'

#major coral groups - tabulate corals, rugose corals, stony corals
tabulata_orders <- c('Favositida', 'Auloporida', 'Heliolitida', 'Sarcinulida', 'Tetradiida', 'Lichenariida', 'Halysitida') #tabulate corals
tabulata_finds <- subset(all_data, all_data$order %in% tabulata_orders)
tabulata_finds$REE_classification <- 'Tabulate_corals'

rugose_orders <- c('Stauriida', 'Cystiphyllida', 'Heterocorallia')
rugose_finds <- subset(all_data, all_data$order %in% rugose_orders)
rugose_finds$REE_classification <- 'Rugose_corals'

scleractinia_order <- c('Scleractinia') #stony corals
scleractinia_finds <- subset(all_data, all_data$order == scleractinia_order)
scleractinia_finds$REE_classification <- 'Stony_corals'

#hydrozoans
hydrozoa_class <- c('Hydrozoa') #hydrozoans
hydrozoa_finds <- subset(all_data, all_data$class == hydrozoa_class)
hydrozoa_finds$REE_classification <- 'Hydrozoans'

#chaetetid sponges 
chaetetids_order <- c('Chaetetida') #sponge
chaetetid_finds <- subset(all_data, all_data$order == chaetetids_order)
chaetetid_finds$REE_classification <- 'Chaetetids'

#reef building tube worms
worm_families <- c('Sabellariidae', 'Serpulidae') #tube worms - honecomb worms and chrittmas tree worms
worm_finds <- subset(all_data, all_data$family %in% worm_families)
worm_finds$REE_classification <- 'Tube_worms'


all_reef_builders <- rbind(archaeocyath_finds,
                           stromatoporoid_finds,
                           hexactinellida_finds,
                           rudist_finds,
                           tabulata_finds,
                           rugose_finds,
                           scleractinia_finds,
                           hydrozoa_finds,
                           chaetetid_finds,
                           worm_finds
                           )

save(all_reef_builders, file='Data/Reef_Ecosystem_Engineers.RData')



