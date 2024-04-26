#Assemble dataset from most recent PBDB download:
#NOTE - WILL NOT REPRODUCE EXACT RESULTS FROM MANUSCRIPT, 
#TO REPRODUCE MANUSCRIPT RESULTS, USE ALREADY EXISTING PHANEROZOIC_CLEAN_FINAL.RDATA IN THE GITHUB REPOSITORY 

#Most of this uses Kocsis et al. 2019 - cite accordingly :) 

library(divDyn)
data(stages)
data(keys)

#=== load pbdb data ===#
filename <- '' #entire filename here 
dat <- as.data.frame(read_csv(filename))
og_dat <- dat #save to keep track of how much data lost at each step  

#=== EARLY PALEOZOIC DATA LOADING - ddPhanero FROM KOCSIS ET AL 2019===#
load(url(
  "https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/cambStrat.RData"))

load(url(
  "https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/ordStrat.RData"))

#=== FUNCTION LOADING ===#
source('Assembly/pbdb_binning_functions.R')
source('Assembly/pbdb_cleaning_functions.R')

#=== Taxonomic cleaning - entire Phanerozoic ====#
#You will lose a large chunk of data here - don't worry! 
dat <- clean_marine(dat) 
clean_dat <- dat #save to keep track of this 

#==== Cambrian stratigraphic binning ====#
source(
  "https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2018-08-31/cambProcess.R")
Cambrian_binned <- dat #saving to keep track
table(dat$stg) #check this worked - should be up to stg 13

#==== Ordovician stratigraphic binning ====#
source(
  "https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2019-05-31/ordProcess.R")
Ordovician_Cambrian_binned <- dat #saving to keep track
table(dat$stg) #check this worked - should be up to stg 27

#==== stratigraphic binning for the rest of Phanerozoic =====#
dat[is.na(dat$late_interval),'late_interval'] <- "" #change empty late intervals to NAs 
dat <- stage_binning(dat)

#let's check that worked -- sometimes because of the way the PBDB may have been download or read in, you may only have ~20k values. This means something has gone wrong! 
table(dat$stg) #check - should be up to stg 95 now
sum(as.data.frame(table(dat$stg))$Freq) #and this should still be a large number -- compare to nrow(Ordovician_Cambrian_binned) to make sure you haven't lost an unreasonable amount of data

#and now we remove any data without a stage assignment
cleaned_binned_dat <- dat #save to compare how much data you are about to lose 
dat <- dat[!is.na(dat$stg),] 
nrow(dat) #this should match what you get for line 49

#and finally put names to stgs
#this can take a while -- your dataframe is still big! 
dat$stage <- rep(NA, nrow(dat))
for(i in 1:length(dat$stg)){
  #find the stage number assinged from early and late interval
  stage_number <- dat$stg[i]
  #what is the corresponding stage in stages?
  stagename <- stages$stage[stage_number]
  #put that back in the dataset
  dat$stage[i] <- stagename
}

#==== formation cleaning - from these authors ====#
formation_cleaning <- read.csv('Assembly/formation_sorting.csv')

colnames(dat)[which(colnames(dat)=='formation')] <- 'old_formation' #preserve original formation assignments
dat <- subset(dat !(old_formation=='')) #remove any data without a formation assignment 

dat$formation <- rep(NA, nrow(dat)) #create new column to save clean formation assignments 
for(i in 1:nrow(dat)){
  this.old.formation <- dat$old_formation[i] #what formation is this occurrence in?
  corrected.formation.row <- subset(formation_cleaning, old_formation==this.old.formation) #look up in formation_cleaning
  corrected.formation <- corrected.formation.row$change_to #save correct formation
  if(length(corrected.formation>0)){ #subtleties in characters and character case that are difficult to catch with the eye can return multiple correct formations, but it is okay to pick the first one from the list (since we are consistent with this, things get groups together as they should)
    dat$formation[i] <- corrected.formation[1]
  }
}

all_data <- dat
save(all_data, file='Data/Phanerozoic_clean_final.RData')

