#Author: Alison Cribb
#Summary: Calculates effect size of bioturbating ecosystem engineers on generic richness, H, and J for each stage through the Phanerozoic

set.seed(541)

#=== USER INPUTS ===#
#how would you like to treat formation subsampling?
#form.subsampling <- 'none'        #no formation subsampling
form.subsampling <- 'occurrences' #subsampling occurrences per formation
#form.subsampling <- 'collections' #subsampling collections per formation

#==== Packages =====#
library(divDyn)

#===== Data input ======#
load('Data/Phanerozoic_clean_final.RData') #phanerozoic PBDB data (date accessed: 1 November 2023)
all_data <- subset(all_data, !(is.na(formation))) #remove data without formation assignments
all_data <- subset(all_data, !(formation=='')) #remove data without formation assignments

load('Data/Bioturbators_data.RData')
ecoeng_genera <- unique(bioturbators_data$genus) #get each ecosystem engineering genus name
ecoeng_formations <- unique(bioturbators_data$formation) #get each formation name containing ecosystem engineers

data(stages) #stage info from divDyn
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]

#=== set up for analysis ===#
#set subsampling variables
n.quota <- 700 #how many fossils per stage?
occs.n.forms <- 20 #how many fossils per formation?
colls.n.forms <- 5 #how many collections per formation?
iter <- 1000 #how many iterations?

#set up dataframe for output
#genrich = generic richness
#H = Shannon's Diversity
#J = Evenness (Pielou Index)
variables <- c('period', 'stage', 'mid_ma', 'n_EE_forms', 'n_nonEE_forms',
               'M1_genrich', 'M1_genrich_sd', 'M2_genrich', 'M2_genrich_sd', 'HedgesG_genrich', 'g_genrich_sd',
               'M1_H', 'M1_H_sd', 'M2_H', 'M2_H_sd', 'HedgesG_H', 'g_H_sd',
               'M1_dom', 'M1_dom_sd', 'M2_dom', 'M2_dom_sd', 'HedgesG_Dominance', 'g_dom_sd')
results_df <- as.data.frame(matrix(NA, 
                                   nrow=length(stage_names),
                                   ncol=length(variables)))
colnames(results_df) <- variables

#=========================#
#====     ANALYSIS    ====#
#=========================#
for(i in 1:length(stage_names)){
  
  #time data
  this.stage <- stage_names[i] 
  this.mid <- stages[which(stages$stage==this.stage),'mid']
  this.period <- stages[which(stages$stage==this.stage),'system']
  results_df$stage[i] <- this.stage
  results_df$mid_ma[i] <- this.mid
  results_df$period[i] <- this.period
  
  #get data for each stage
  this.stage.data <- subset(all_data, stage==this.stage)

  #ecosystem engineering data
  all_EE_data <- subset(this.stage.data, formation %in% ecoeng_formations) #presence data for entire stage
  results_df$n_EE_forms[i] <- length(unique(all_EE_data$formation)) #how many formations have ecosystem engineers in this stage?
  all_nonEE_data <- subset(this.stage.data, !(formation %in% ecoeng_formations)) #absence data for entire stage
  results_df$n_nonEE_forms[i] <- length(unique(all_nonEE_data$formation)) #how many formations do not have ecosystem engineers in this stage?
  
  #temporary vectors to save subsampled data 
  M1_genrich_iters <- rep(NA, iter)
  M2_genrich_iters <- rep(NA, iter)
  HedgesG_genrich_iters <- rep(NA, iter)
  
  M1_H_iters <- rep(NA, iter)
  M2_H_iters <- rep(NA, iter)
  HedgesG_H_iters <- rep(NA, iter)
  
  M1_dom_iters <- rep(NA, iter)
  M2_dom_iters <- rep(NA, iter)
  HedgesG_dom_iters <- rep(NA, iter)
  
  #subsampling
  for(j in 1:iter){
    
    #subsample/bootstrap n.quota occurrences per stage
    max.samples <- nrow(this.stage.data)
    row_idxs <- sample(max.samples, n.quota, replace=TRUE)
    subbed.data <- this.stage.data[row_idxs,]
    
    #present ecosystem engineering data from stage-subsampled occurrences:
    presence_data <- subset(subbed.data, formation %in% ecoeng_formations)
    
    #**** If you only want to consider ecosystem engineer impacts on non-ecosystem engineering taxa, but see manuscript for discussion for why we do not do this in the final analyses, change to:
    #presence_data_all <- subset(subbed.data, formation %in% ecoeng_formations)
    #presence_data <- subset(presence_data_all, !(genus %in% ecoeng_genera))
    
    #absent ecosystem engineering data from stage-subsampled occurrences:
    absence_data <- subset(subbed.data, !(formation %in% ecoeng_formations))
    
    #get lists of presence and absence formations from that subsampled data 
    presence_formations <- unique(presence_data$formation)
    absence_formations <- unique(absence_data$formation)
    
    #=== Collect ecological statistics: ===#
    #n1,2 = sample size (no. formations)
    #x1,2 = mean generic richness/diversity/evenness (between formations)
    #s1,2 = standard deviation of means 
    
    n1 <- length(presence_formations) #number of formations containing bioturbating ecosystem engineers
    n2 <- length(absence_formations) #number of formations NOT containing bioturbating ecosystem engineers
    
    #== presence statistics (x1) ==#
    genrich.presence.temp <- rep(NA, n1)
    H.presence.temp <- rep(NA, n1) #set up temporary vector to save Shannon's Diversity for each of the n1 presence formations 
    Dom.presence.temp <- rep(NA, n1) #and for Simpson's dominance 
    for(k in 1:n1){
      this.formation <- presence_formations[k] 
      this.formation.data <- subset(presence_data, formation==this.formation) #get data for each formation
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling occurrences per formation,
        form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE)  #get random rows
        this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }

      if(form.subsampling=='collections'){ #if subsampling collections per formation
        form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
        sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get random numbers that correspond to places in form_colls
        sub.colls <- form_colls[sub.idxs] #and pull those random collection numbers out of form_colls
        this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those collections 
      }
      
      this.presence_genera <- unique(this.formation.subbed$genus) #get list of all genera in the subsampled formation
      this.presence_abundance_data <- as.data.frame(matrix(NA, 
                                                           nrow=length(this.presence_genera),
                                                           ncol=2)) #set up to collect how many of each genus there are 
      colnames(this.presence_abundance_data) <- c('gen', 'n')
      for(l in 1:nrow(this.presence_abundance_data)){ #for each genus
        this.genus <- this.presence_genera[l] 
        this.presence_abundance_data$gen[l] <- this.genus
        this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
        this.presence_abundance_data$n[l] <- nrow(this.genus.data) #count how many occurrences there are 
      }
      genrich.presence.temp[k] <- nrow(this.presence_abundance_data) #calculate generic richness from n. taxa (n. rows) from abundance matrix
      
      #Shannon's diversity (H)
      presence_tot <- sum(this.presence_abundance_data$n) #total number of occurrences of all genera (formation size)
      presence_gen_props <- this.presence_abundance_data$n/presence_tot #relative abundance for each genus
      presence_shannon_div <- -sum(presence_gen_props*log(presence_gen_props)) #Shannon's Diversity formula 
      H.presence.temp[k] <- presence_shannon_div #and save
      # was this originally to add points to an ongoing plot? -WG
      #points(nrow(this.presence_abundance_data), presence_shannon_div)
      
      #and we can use this all to calculate Simpson's dominance
      D.presence <- (sum(presence_gen_props^2)) #sum of squared generic proportions 
      Dom.presence.temp[k] <- 1/D.presence #Simpson's dominance=1/D
      
    }
    
    #Generic richness presence statistics
    x1 <- mean(genrich.presence.temp, na.rm=TRUE) #find mean generic richness across all of the presence formations
    s1 <- sd(genrich.presence.temp, na.rm=TRUE) #and the standard deviation
    #save
    M1_genrich_iters[j] <- x1
    
    #Shannon's Diversity presence statistics 
    x1_H <- mean(H.presence.temp, na.rm=TRUE) #mean Shannon's Diversity across all presence formations
    s1_H <- sd(H.presence.temp, na.rm=TRUE) #and the standard deviation
    M1_H_iters[j] <- x1_H #save
    
    #Simpson's Dominance statistics 
    x1_dom <- mean(Dom.presence.temp, na.rm=TRUE)
    s1_dom <- sd(Dom.presence.temp, na.rm=TRUE)
    M1_dom_iters[j] <- x1_dom #save 
    
    
    #== absence statistics (x2) ==#
    genrich.absence.temp <- rep(NA, n2)
    H.absence.temp <- rep(NA, n2) #set up temporary vector to save Shannon's Diversity for each of the n2 absence formations 
    Dom.absence.temp <- rep(NA, n2) #and for Simpson's dominance 
    for(k in 1:n2){
      this.formation <- absence_formations[k]
      this.formation.data <- subset(absence_data, formation==this.formation) #get data for each formation
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling occurrences per formation,
        form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE) #get random rows
        this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }
       
      if(form.subsampling=='collections'){  #if subsampling collections per formation
        form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
        sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get random numbers that correspond to places in form_colls
        sub.colls <- form_colls[sub.idxs] #and pull those random collection numbers out of form_colls
        this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those collections 
      }
      
      this.absence_genera <- unique(this.formation.subbed$genus) #get list of all of the genera in this absence formation
      this.absence_abundance_data <- as.data.frame(matrix(NA, 
                                                          nrow=length(this.absence_genera),
                                                          ncol=2)) #set up to collection how many occurrences of each genus in the formation
      colnames(this.absence_abundance_data) <- c('gen', 'n') 
      for(l in 1:nrow(this.absence_abundance_data)){ #for each genus
        this.genus <- this.absence_genera[l] 
        this.absence_abundance_data$gen[l] <- this.genus
        this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
        this.absence_abundance_data$n[l] <- nrow(this.genus.data)  #count how many occurrences there are 
      }
      genrich.absence.temp[k] <- nrow(this.absence_abundance_data) #calculate generic richness from n. taxa (n. rows) from abundance matrix 
      
      #Shannon's diversity 
      absence_tot <- sum(this.absence_abundance_data$n) #total number of occurrences of all genera (formation size)
      absence_gen_props <- this.absence_abundance_data$n/absence_tot #relative abundance for each genus
      absence_shannon_div <- -sum(absence_gen_props*log(absence_gen_props)) #Shannon's Diversity Formula
      H.absence.temp[k] <- absence_shannon_div
      
      #and we can sue this all to calculate Simpson's dominance
      D.absence <- (sum(absence_gen_props^2)) #sum of squared generic proportions 
      Dom.absence.temp[k] <- 1/D.absence #Simpson's dominance=1/D
    }
    

    #Generic richness absence statistics
    x2 <- mean(genrich.absence.temp, na.rm=TRUE) #find mean generic richness across all of the absence formations
    s2 <- sd(genrich.absence.temp, na.rm=TRUE) #and the standard deviation 
    M2_genrich_iters[j] <- x2 #save

    #Shannon's Diversity absence statistics
    x2_H <- mean(H.absence.temp, na.rm=TRUE) #mean Shannon's Diversity across all absence formations
    s2_H <- sd(H.absence.temp, na.rm=TRUE) #and the standard deviation
    M2_H_iters[j] <- x2_H #save
    
    #Simpson's Dominance absence statistics 
    x2_dom <- mean(Dom.absence.temp, na.rm=TRUE)
    s2_dom <- sd(Dom.absence.temp, na.rm=TRUE)
    M2_dom_iters[j] <- x2_dom #save 
    
    #== Effect sizes ==#
    #Effect size for generic richness
    genrich_g <- (x1-x2)/( sqrt( ( ((n1-1)*(s1^2)) + ((n2-1)*(s2^2)) ) / (n1+n2-2)   )  ) #calculate Hedges G comparing the presence and absence data
    HedgesG_genrich_iters[j] <- genrich_g #save 
 
    #Effect size for Shannon's Diversity 
    shannondiv_g <- (x1_H-x2_H)/( sqrt( ( ((n1-1)*(s1_H^2)) + ((n2-1)*(s2_H^2)) ) / (n1+n2-2)   )  ) #Hedges G comparing Shannon's Diversity of presence and absence data
    HedgesG_H_iters[j] <- shannondiv_g #save 

    #Effect size for Simpson's Dominance 
    dominance_g <- (x1_dom-x2_dom)/( sqrt( ( ((n1-1)*(s1_dom^2)) + ((n2-1)*(s2_dom^2)) ) / (n1+n2-2)   )  )
    HedgesG_dom_iters[j] <- dominance_g

  }
  
  #save data with errors in main results 
  #Generic richness
  # it might make more sense to do weighted means and weighted standard deviations here
  # this will take into account the standard deviations associated with each mean and will downweight those with larger errors
  # of course you can't do this for the variables where you don't have a SD/SE for each iteration (e.g., HedgesG)
  # see Hmisc::wtd.mean and Hmisc::wtd.var (you want the sqrt of the latter)
  # I usually use 1/(SE^2) for the weights (and you need to use normwt=TRUE)
  # -WG
  results_df$M1_genrich[i] <- mean(M1_genrich_iters, na.rm=TRUE)
  results_df$M1_genrich_sd[i] <- sd(M1_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich[i] <- mean(M2_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich_sd[i] <- sd(M2_genrich_iters, na.rm=TRUE)
  results_df$HedgesG_genrich[i] <- mean(HedgesG_genrich_iters, na.rm=TRUE)
  results_df$g_genrich_sd[i] <- sd(HedgesG_genrich_iters, na.rm=TRUE)
  
  #Shannon's diversity
  results_df$M1_H[i] <- mean(M1_H_iters, na.rm=TRUE)
  results_df$M1_H_sd[i] <- sd(M1_H_iters, na.rm=TRUE)
  results_df$M2_H[i] <- mean(M2_H_iters, na.rm=TRUE)
  results_df$M2_H_sd[i] <- sd(M2_H_iters, na.rm=TRUE)
  results_df$HedgesG_H[i] <- mean(HedgesG_H_iters, na.rm=TRUE)
  results_df$g_H_sd[i] <- sd(HedgesG_H_iters, na.rm=TRUE)
  
  #Dominance 
  results_df$M1_dom[i] <- mean(M1_dom_iters, na.rm=TRUE)
  results_df$M1_dom_sd[i] <- sd(M1_dom_iters, na.rm=TRUE)
  results_df$M2_dom[i] <- mean(M2_dom_iters, na.rm=TRUE)
  results_df$M2_dom_sd[i] <- sd(M2_dom_iters, na.rm=TRUE)
  results_df$HedgesG_Dominance[i] <- mean(HedgesG_dom_iters, na.rm=TRUE)
  results_df$g_dom_sd[i] <- sd(HedgesG_dom_iters, na.rm=TRUE)
  
  print(paste('finished stage:', stage_names[i]))
  
}

#Save results 
bioturbation_results_df <- results_df

# doesn't work if the Output folder doesn't exist yet -WG
#save(bioturbation_results_df, file='Output/effectsizes_bioturbation_collsub.RData')
save(bioturbation_results_df, file='Output/effectsizes_bioturbation_occsub.RData')
#save(bioturbation_results_df, file='Output/effectsizes_bioturbation_noformsub.RData')


