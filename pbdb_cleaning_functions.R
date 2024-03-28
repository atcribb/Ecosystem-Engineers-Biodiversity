#Functions to clean marine pbdb data
#Author: Alison Cribb (A.T.Cribb@soton.ac.uk)
#Created: 20 September 2023
#Last edited: 1 November 2023


#======== PBDB URL BOX DO NOT TOUCH =========#
#Marine:
#data1.2/occs/list.csv?idqual=certain&pres=regular&interval=Cambrian,Holocene&envtype=marine&pgm=gplates,scotese,seton&show=full,img,lith,env
#=============================================#

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










