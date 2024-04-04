#Author: Alison Cribb
#Summary: Calculates effect size of bioturbating ecosystem engineers on generic richness, H, and J for each stage through the Phanerozoic

set.seed(541)

#=== USER INPUTS ===#
#how would you like to treat formation subsampling?
#form.subsampling <- 'none'        #no formation subsampling
form.subsampling <- 'occurrences' #subsampling - 20 occurrences per formation
#form.subsampling <- 'collections' #subsampling - 5 collections per formation

#==== Packages =====#
library(divDyn)

#===== Data input ======#
load('Phanerozoic_clean_final.RData') #phanerozoic PBDB data (date accessed: 1 November 2023)
all_data <- subset(all_data, !(is.na(formation))) #remove data without formation assignments
all_data <- subset(all_data, !(formation=='')) #remove data without formation assignments

load('Bioturbators_data.RData')
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
#genrich = generich richness
#H = Shannon's Diversity
#J = Evenness (Pielou Index)
variables <- c('period', 'stage', 'mid_ma', 'n_EE_forms', 'n_nonEE_forms',
               'M1_genrich', 'M1_genrich_sd', 'M2_genrich', 'M2_genrich_sd', 'HedgesG_genrich', 'g_genrich_sd',
               'M1_H', 'M1_H_sd', 'M2_H', 'M2_H_sd', 'HedgesG_H', 'g_H_sd',
               'M1_dom', 'M1_dom_sd', 'M2_dom', 'M2_dom_sd', 'HedgesG_Dominance', 'g_dom_sd',
               'M1_J', 'M1_J_sd', 'M2_J', 'M2_J_sd', 'HedgesG_J', 'g_J_sd')
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
  
  # M1_J_iters <- rep(NA, iter)
  # M2_J_iters <- rep(NA, iter)
  # HedgesG_J_iters <- rep(NA, iter)
  
  M1_dom_iters <- rep(NA, iter)
  M2_dom_iters <- rep(NA, iter)
  HedgesG_dom_iters <- rep(NA, iter)
  
  #subsampling
  for(j in 1:iter){
    
    #subsample n.quota occurrences per stage 
    max.samples <- nrow(this.stage.data)
    row_idxs <- sample(max.samples, n.quota, replace=TRUE)
    subbed.data <- this.stage.data[row_idxs,]
    
    #present ecosystem engineering data from stage-subsampled 700 occurrences:
    presence_data <- subset(subbed.data, formation %in% ecoeng_formations)
    
    #**** If you only want to consider ecosystem engineer impacts on non-ecosystem engineering taxa, but see manuscript for discussion for why we do not do this in the final analyses, change to:
    #presence_data_all <- subset(subbed.data, formation %in% ecoeng_formations)
    #presence_data <- subset(presence_data_all, !(genus %in% ecoeng_genera))
    
    #absent ecosystem engineering data from stage-subsampled 700 occurrences:
    absence_data <- subset(subbed.data, !(formation %in% ecoeng_formations))
    
    #get lists of presence and absence formations from that subsampled data 
    presence_formations <- unique(presence_data$formation)
    absence_formations <- unique(absence_data$formation)
    
    #=== Collect ecological statistics: ===#
    #Hedges' G
    #n1,2 = sample size (no. formations)
    #x1,2 = mean generic richness/diversity/evenness (between formations)
    #s1,2 = standard deviation of means 
    
    n1 <- length(unique(presence_formations))
    n2 <- length(unique(absence_formations))
    
    #== GENERIC RICHNESS ==#
    #presence statistics
    genrich.presence.temp <- rep(NA, n1)
    for(k in 1:n1){ #collected generic richness in each n1 formations
      this.formation <- presence_formations[k] #from the list of presence_formations,
      this.formation.data <- subset(presence_data, formation==this.formation) #pick out a single formation
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling 20 occurrences per formation,
      form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE) #get 20 random rows
      this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }
      
      if(form.subsampling=='collections'){ #if subsampling 5 collections per formation
      form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
      sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get 5 random numbers that correspond to palces in form_colls
      sub.colls <- form_colls[sub.idxs] #and pull those 5 random collection numbers out of form_colls
      this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those 5 collections 
      }
      
      genrich.presence.temp[k] <- length(unique(this.formation.subbed$genus)) #calculate generic richness 
    }
    
    x1 <- mean(genrich.presence.temp, na.rm=TRUE) #find mean generic richenss across all of the presence formations
    s1 <- sd(genrich.presence.temp, na.rm=TRUE) #and the standard deviation
    
    #absence statistics
    genrich.absence.temp <- rep(NA, n2)
    for(k in 1:n1){
      this.formation <- absence_formations[k] 
      this.formation.data <- subset(absence_data, formation==this.formation)
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling 20 occurrences per formation,
      form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE) #get 20 random rows
      this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }
      
      if(form.subsampling=='collections'){ #if subsampling 5 collections per formation
      form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
      sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get 5 random numbers that correspond to palces in form_colls
      sub.colls <- form_colls[sub.idxs] #and pull those 5 random collection numbers out of form_colls
      this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those 5 collections 
      }
      
      genrich.absence.temp[k] <- length(unique(this.formation.subbed$genus)) #calculate generic richness
    }
    
    x2 <- mean(genrich.absence.temp, na.rm=TRUE) #find mean generic richness across all of the absence formations
    s2 <- sd(genrich.absence.temp, na.rm=TRUE) #and the standard deviation 
    
    #calculate hedges' g
    genrich_g <- (x1-x2)/( sqrt( ( ((n1-1)*(s1^2)) + ((n2-1)*(s2^2)) ) / (n1+n2-2)   )  ) #calcualte Hedges G comparing the presence and absence data
    
    #save
    M1_genrich_iters[j] <- x1
    M2_genrich_iters[j] <- x2
    HedgesG_genrich_iters[j] <- genrich_g
    
    
    #== SHANNON'S DIVERSITY (H) AND EVENNESS (J) ==#
    #presence statistics (x1)
    H.presence.temp <- rep(NA, n1) #set up temporary vector to save Shannon's Diversit for each of the n1 presence formations 
    Dom.presence.temp <- rep(NA, n1) #and for Simpson's dominance 
    for(k in 1:n1){
      this.formation <- presence_formations[k] 
      this.formation.data <- subset(presence_data, formation==this.formation) #get data for each formation
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling 20 occurrences per formation,
        form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE)  #get 20 random rows
        this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }

      if(form.subsampling=='collections'){ #if subsampling 5 collections per formation
        form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
        sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get 5 random numbers that correspond to places in form_colls
        sub.colls <- form_colls[sub.idxs] #and pull those 5 random collection numbers out of form_colls
        this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those 5 collections 
      }
      
      this.presence_genera <- unique(this.formation.subbed$genus) #get list of all genera in the subsampled formation
      this.presence_abundance_data <- as.data.frame(matrix(NA, 
                                                           nrow=length(this.presence_genera),
                                                           ncol=2)) #set up to collect how many of each genera there are 
      colnames(this.presence_abundance_data) <- c('gen', 'n')
      for(l in 1:nrow(this.presence_abundance_data)){ #for each genus
        this.genus <- this.presence_genera[l] 
        this.presence_abundance_data$gen[l] <- this.genus
        this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
        this.presence_abundance_data$n[l] <- nrow(this.genus.data) #count how many occurrences there are 
      }
      
      #Shannon's diversity (H)
      presence_tot <- sum(this.presence_abundance_data$n) #total number of occurrences of all genera (formation size)
      presence_gen_props <- this.presence_abundance_data$n/presence_tot #relative abundance for each genus
      presence_shannon_div <- -sum(presence_gen_props*log(presence_gen_props)) #Shannon's Diversity formula 
      H.presence.temp[k] <- presence_shannon_div #and save 
      
      #and we can use this all to calculate Simpson's dominance
      D.presence <- (sum(presence_gen_props^2)) #sum of squared generic proprtions 
      Dom.presence.temp[k] <- 1/D.presence #Simpson's dominance=1/D
      
    }
    
    #Shannon's Diversity statistics 
    x1_H <- mean(H.presence.temp, na.rm=TRUE) #mean Shannon's Diversity across all presence formations
    s1_H <- sd(H.presence.temp, na.rm=TRUE) #and the standard deviation
    
    #Simpson's Dominance statistics 
    x1_dom <- mean(Dom.presence.temp, na.rm=TRUE)
    s1_dom <- sd(Dom.presence.temp, na.rm=TRUE)
    
    #absence statistics (x2)
    H.absence.temp <- rep(NA, n2) #set up temporary vector to save Shannon's Diversit for each of the n2 absence formations 
    Dom.absence.temp <- rep(NA, n2) #and for Simpson's dominance 
    for(k in 1:n2){
      this.formation <- absence_formations[k]
      this.formation.data <- subset(absence_data, formation==this.formation) #get data for each formation
      
      if(form.subsampling=='none'){ #if no within formation subsampling, do nothing 
        this.formation.subbed <- this.formation.data
      }
      
      if(form.subsampling=='occurrences'){ #if subsampling 20 occurrences per formation,
        form_row_idxs <- sample(nrow(this.formation.data), occs.n.forms, replace=TRUE) #get 20 random rows
        this.formation.subbed <- this.formation.data[form_row_idxs,] #and pull out that data from this formation data
      }
       
      if(form.subsampling=='collections'){  #if subsampling 5 collections per formation
        form_colls <- unique(this.formation.data$collection_no) #get all of the collection numbers for this formation
        sub.idxs <- sample(length(form_colls), colls.n.forms, replace=TRUE) #get 5 random numbers that correspond to places in form_colls
        sub.colls <- form_colls[sub.idxs] #and pull those 5 random collection numbers out of form_colls
        this.formation.subbed <- subset(this.formation.data, collection_no %in% sub.colls) #and then pull all of the data for those 5 collections 
      }
      
      this.absence_genera <- unique(this.formation.subbed$genus) #get list of all of the genera in this absence formation
      this.absence_abundance_data <- as.data.frame(matrix(NA, 
                                                          nrow=length(this.absence_genera),
                                                          ncol=2)) #set up to collection how many occurrences of each genera in the formation
      colnames(this.absence_abundance_data) <- c('gen', 'n') 
      for(l in 1:nrow(this.absence_abundance_data)){ #for each genus
        this.genus <- this.absence_genera[l] 
        this.absence_abundance_data$gen[l] <- this.genus
        this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
        this.absence_abundance_data$n[l] <- nrow(this.genus.data)  #count how many occurrences there are 
      }
      
      #Shannon's diversity 
      absence_tot <- sum(this.absence_abundance_data$n) #total number of occurrences of all genera (formation size)
      absence_gen_props <- this.absence_abundance_data$n/absence_tot #relative abundance for each genus
      absence_shannon_div <- -sum(absence_gen_props*log(absence_gen_props)) #Shannon's Diversity Formula
      H.absence.temp[k] <- absence_shannon_div
      
      #and we can sue this all to calculate Simpson's dominance
      D.absence <- (sum(absence_gen_props^2)) #sum of squared generic proprtions 
      Dom.absence.temp[k] <- 1/D.absence #Simpson's dominance=1/D
      
    }
    
    #Shannon's Diversity presence statistics
    x2_H <- mean(H.absence.temp, na.rm=TRUE) #mean Shannon's Diversity across all absence formations
    s2_H <- sd(H.absence.temp, na.rm=TRUE) #and the standard deviation
    
    #Effect size for Shannon's Diversity 
    shannondiv_g <- (x1_H-x2_H)/( sqrt( ( ((n1-1)*(s1_H^2)) + ((n2-1)*(s2_H^2)) ) / (n1+n2-2)   )  ) #Hedges G comparing Shannon's Diversity of presence and absence data

    #Save 
    M1_H_iters[j] <- x1_H
    M2_H_iters[j] <- x2_H
    HedgesG_H_iters[j] <- shannondiv_g
    
    #Simpson's Dominance absence statistics 
    x2_dom <- mean(Dom.absence.temp, na.rm=TRUE)
    s2_dom <- sd(Dom.absence.temp, na.rm=TRUE)
    
    #Effect size for Simpson's Dominance 
    dominance_g <- (x1_dom-x2_dom)/( sqrt( ( ((n1-1)*(s1_dom^2)) + ((n2-1)*(s2_dom^2)) ) / (n1+n2-2)   )  )
    
    #Save
    M1_dom_iters[j] <- x1_dom
    M2_dom_iters[j] <- x2_dom
    HedgesG_dom_iters[j] <- dominance_g
    
    
    # #==== EVENNESS (J) ====#
    # #Note that for this metric, we are really analyzing how common very diverse formations are 
    #     #When J is closer to 1, more formations are closer to Hmax
    #     #When J is closer to 0, the only a few formations are close to Hmax 
    # #J = H_i/max(H)
    # #presence statistics
    # max_H_presence <- max(H.presence.temp) #get maximum Shannon's Diversity index across presence formations 
    # J.presence.temp <- H.presence.temp/max_H_presence #get evenness indices (J) for each formation
    # 
    # x1_J <- mean(J.presence.temp, na.rm=TRUE) #get mean J across the formations 
    # s1_J <- sd(J.presence.temp, na.rm=TRUE) #and standard deviation 
    # 
    # #absence statistics 
    # max_H_absence <- max(H.absence.temp) #get maximum Shannon's Diversity index across asbence formations 
    # J.absence.temp <- H.absence.temp/max_H_absence #get evenness indices (J) for each formation
    # 
    # x2_J <- mean(J.absence.temp, na.rm=TRUE) #get mean J across the formations 
    # s2_J <- sd(J.absence.temp, na.rm=TRUE) #and standard deviation 
    # 
    # evenness_g <- (x1_J-x2_J)/( sqrt( ( ((n1-1)*(s1_J^2)) + ((n2-1)*(s2_J^2)) ) / (n1+n2-2)   )  ) #and calculate Hedges' g comparing Pielout's index in presence and absence data
    
    # #save
    # M1_J_iters[j] <- x1_J
    # M2_J_iters[j] <- x2_J
    # HedgesG_J_iters[j] <- evenness_g
    # 
  }
  
  #print('test point 7')
  
  #save data with errors in main results 
  #Generic richness
  results_df$M1_genrich[i] <- mean(M1_genrich_iters, na.rm=TRUE)
  results_df$M1_genrich_sd[i] <- sd(M1_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich[i] <- mean(M2_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich_sd[i] <- sd(M2_genrich_iters, na.rm=TRUE)
  results_df$HedgesG_genrich[i] <- mean(HedgesG_genrich_iters, na.rm=TRUE)
  results_df$g_genrich_sd[i] <- sd(HedgesG_genrich_iters, na.rm=TRUE)
  
  #Shannon's diversity
  results_df$M1_H[i] <- mean(M1_H_iters, na.rm=TRUE)
  results_df$M1_H_sd[i] <- sd(M1_H_iters,)
  results_df$M2_H[i] <- mean(M2_H_iters, na.rm=TRUE)
  results_df$M2_H_sd[i] <- sd(M2_H_iters, na.rm=TRUE)
  results_df$HedgesG_H[i] <- mean(HedgesG_H_iters, na.rm=TRUE)
  results_df$g_H_sd[i] <- sd(HedgesG_H_iters, na.rm=TRUE)
  
  # #Evenness
  # results_df$M1_J[i] <- mean(M1_J_iters, na.rm=TRUE)
  # results_df$M1_J_sd[i] <- sd(M1_J_iters, na.rm=TRUE)
  # results_df$M2_J[i] <- mean(M2_J_iters, na.rm=TRUE)
  # results_df$M2_J_sd[i] <- sd(M2_J_iters, na.rm=TRUE)
  # results_df$HedgesG_J[i] <- mean(HedgesG_J_iters, na.rm=TRUE)
  # results_df$g_J_sd[i] <- sd(HedgesG_J_iters, na.rm=TRUE)
  
  #Dominance 
  results_df$M1_dom[i] <- mean(M1_dom_iters, na.rm=TRUE)
  results_df$M1_dom_sd[i] <- sd(M1_dom_iters, na.rm=TRUE)
  results_df$M2_dom[i] <- mean(M2_dom_iters, na.rm=TRUE)
  results_df$M2_dom_sd[i] <- sd(M2_dom_iters, na.rm=TRUE)
  results_df$HedgesG_Dominance[i] <- mean(HedgesG_dom_iters, na.rm=TRUE)
  results_df$g_dom_sd[i] <- sd(HedgesG_dom_iters, na.rm=TRUE)
  
  
  print(paste('finished stage:', stage_names[i]))
  
}

bioturbation_results_df <- results_df

save(bioturbation_results_df, file='effectsizes_bioturbation.RData')


