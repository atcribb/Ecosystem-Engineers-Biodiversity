#Functions to clean marine and terrestrial pbdb data

#FOLLOWING KOCSIS ET AL 2019 -- CITE THEM ACCORDINGLY! :)

library(divDyn)
library(dplyr)

clean_marine <- function(data){
  
  
  #Following ddPhanero (Kocscis et al.)
  #==== Taxonomic Filtering =====#
  #omit occurrences not identified to genus level
  data <- subset(data, accepted_rank %in% c('genus', 'species'))
  data <- data[data$genus!='',] 

  #phyla-level filtering
  marine.phyla <- c("",
                     "Agmata",
                     "Annelida",
                     "Bilateralomorpha",
                     "Brachiopoda",
                     "Bryozoa",
                     "Calcispongea",
                     "Chaetognatha",
                     "Cnidaria",
                     "Ctenophora",
                     "Echinodermata",
                     "Entoprocta",
                     "Foraminifera",
                     "Hemichordata",
                     "Hyolitha",
                     "Mollusca",
                     "Nematoda",
                     "Nematomorpha",
                     "Nemertina",
                     "Onychophora",
                     "Petalonamae",
                     "Phoronida",
                     "Platyhelminthes",
                     "Porifera",
                     "Rhizopodea",
                     "Rotifera",
                     "Sarcomastigophora",
                     "Sipuncula",
                     "Uncertain",
                     "Vetulicolia",
                     "")
  select.phyla <- (data$phylum %in% marine.phyla) #logical vector of where phyla of interest are
  #data <- data[select.phyla,]

  #class-level filtering
  marine.class <- c(
                    "Acanthodii",
                    "Actinopteri",
                    "Actinopterygii",
                    "Agnatha",
                    "Cephalaspidomorphi",
                    "Chondrichthyes",
                    "Cladistia",
                    "Coelacanthimorpha",
                    "Conodonta",
                    "Galeaspida",
                    "Myxini",
                    "Osteichthyes",
                    "Petromyzontida",
                    "Plagiostomi",
                    "Pteraspidomorphi",
                    "Artiopoda",
                    "Branchiopoda",
                    "Cephalocarida",
                    "Copepoda",
                    "Malacostraca",
                    "Maxillopoda",
                    "Megacheira",
                    "Merostomoidea",
                    "Ostracoda",
                    "Paratrilobita",
                    "Pycnogonida",
                    "Remipedia",
                    "Thylacocephala",
                    "Trilobita",
                    "Xiphosura"
                  )
  select.class <- (data$class %in% marine.class) #logical vector of where classes of interest are 
  #data <- data[select.class,]

  #filtering for mammals
  marine.mammals.order <- c('Cetacea', 'Sirenia')
  select.mammals.order <- (data$order %in% marine.mammals.order)
  
  #and mammalian carnivores 
  marine.mammals.family <- c("Otariidae", "Phocidae", "Desmatophocidae")
  select.mammals.family <- (data$family %in% marine.mammals.family)
  
  #filtering for marine reptiles
  marine.reptiles.order <- c("Eosauropterygia",
                               "Hupehsuchia",
                               "Ichthyosauria",
                               "Placodontia",
                               "Sauropterygia",
                               "Thalattosauria")
  select.marine.reptiles <- (data$order %in% marine.reptiles.order)
  
  #filtering for sea turtles
  turtles.family <- c("Cheloniidae",
                       "Protostegidae",
                       "Dermochelyidae",
                       "Dermochelyoidae",
                       "Toxochelyidae",
                       "Pancheloniidae")
  select.turtles <- (data$family %in% turtles.family)
  
  #finally, subset this data with the multiple filters
  data <- data[select.phyla | select.class | select.mammals.order | select.mammals.family | select.marine.reptiles | select.turtles ,]

  #resolve homonymies by combining class names and genus names to create individual entries
  data$clgen <- paste(data$class, data$genus)
         
  

  #==== Environmental Filtering =====#     
  #remove occurrences in environments that are more likely to be terrestrial taxa 
  omit.environments <- c(
    "\"floodplain\"", "alluvial fan", "cave", "\"channel\"", "channel lag" , 
    "coarse channel fill", "crater lake", "crevasse splay", "dry floodplain", 
    "delta plain", "dune", "eolian indet.", "fine channel fill", "fissure fill", 
    "fluvial indet.", "fluvial-lacustrine indet.", "fluvial-deltaic indet.", 
    "glacial", "interdune", "karst indet.", "lacustrine - large", 
    "lacustrine - small", "lacustrine delta front", "lacustrine delta plain", 
    "lacustrine deltaic indet.", "lacustrine indet.", 
    "lacustrine interdistributary bay", "lacustrine prodelta", "levee", "loess", 
    "mire/swamp", "pond", "sinkhole", "spring", "tar", "terrestrial indet.", 
    "wet floodplain")
  data <- subset(data, !(environment %in% omit.environments)) #removal

  return(data)

  
}

# #==== test ====#
# load('test_data.RData')
# marine_data <- clean_marine(test.data)




#======== TERRESTRIAL CLEANING FUNCTION ========#
####===========================  CHECK YOUR DATA BEFORE YOU BEGIN  ===========================######
#     For things to go smoothly, when you downloaded your PBDB data did you....                    #
#          - select "terrestrial" in Select by geologic context>Environment ?                      #
#          - select "regular taxa only" in Select by taxonomy>Preservation ?                       #
#          - deselect "include metadata at the beginning of the output" in Choose output options ? #
#==================================================================================================#

#make sure trace_removals and egg_removals .RData files are loaded (these are from E. Dunne)
load('trace_removal_terms.RData')
load('egg_removal_terms.RData')
load('marine_removal_terms.RData')

#get marine genera from your clean marine dataset
#marine_genus <- unique(marine_data$genus)

clean_terrestrial <- function(data){
  
  #Following Cribb, Formoso et al. (2023) and Emma Dunne (github.com/emmadunne/pbdb_clean)
  
  #==== Taxonomic Filtering =====#
  #omit occurrences not identified to genus level
  data <- subset(data, accepted_rank %in% c('genus', 'species'))
  data <- data[data$genus!='',]
  
  #omit plants, insects, microbes, and meiofauna
  terrestrial.phyla <- c('Chordata') #for future, add in Arthropoda and Mollusca for insects and inverts
  data <- subset(data, phylum %in% terrestrial.phyla)
  
  omit.fishes <- c('Actinopteri', 'Chondrichthyes', 'Osteichthyes')
  data <- subset(data, !(class %in% omit.fishes))
  
  #omit occurrences attributed to trace fossils, eggs, marine taxa, and wastebasket taxa (from Emma Dunne)
  wastebasket_removals <- c('Crocodylus', 'Alligator', 'Testudo', 'Lacerta', 'Fenestrosaurus', 'Ovoraptor', 'Ornithoides', 'Placodermi') #remove wastebasket taxa to improve taxonomic resolution
  removal.terms <- c(trace_removals, egg_removals, wastebasket_removals) #marine_removals
  data <- subset(data, !(class %in% removal.terms))
  data <- subset(data, !(order %in% removal.terms))
  data <- subset(data, !(family %in% removal.terms))
  data <- subset(data, !(genus %in% removal.terms))
  
  #omit marine genear that end up in the terrestrial dataset
  #data <- subset(data, !(genus %in% marine_genus))
  
  #omit any lingering trace fossils 
  data <- subset(data, !(pres_mode=='trace'))
  
  #omit iffy taxonomic assignments to resolve taxonomic resolution
  data <- data %>% filter(!grepl("cf\\.|aff\\.|\\?|ex\\. gr\\.|sensu lato|informal|\\\"", identified_name))
  
  #===== Environmental filtering =====#
  omit.environments <- c('basinal (carbonate)',             #remove environments least likely to be impacted by bloat and float 
                         'basinal (siliceous)',
                         'basinal (siliciclastic)',
                         'carbonate indet.',
                         'deep-water indet.', 
                         'marine indet.',
                         'open shallow subtidal', 
                         'perireef or subreef', 'peritidal', 'reef, buildup or bioherm', 
                         'shoreface', 'slope', 'slope/ramp reef', 'transition zone/lower'
                         )  
  data <- subset(data, !(environment %in% omit.environments))
  
  return(data)
  
}











