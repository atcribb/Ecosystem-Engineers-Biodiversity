#Authors: Alison Cribb, Will Gearty
#Summary: Assess sampling biases influence on bioturbation effect size results

#=========#
# clear old data
rm(list = ls())

library(divDyn)
data("stages", package="divDyn")
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]
period_names <- unique(stages[which(stages$stage %in% stage_names), 'system'])
period.cols <- unique(stages[which(stages$stage %in% stage_names), 'systemCol'])
library(ggplot2)
library(deeptime)

sampling_theme <- theme(
  legend.position="inside",
  legend.position.inside=c(0.85, 0.9),
  legend.title=element_blank(),
  plot.title=element_text(hjust=0.5),
  axis.text = element_text(color = "black"),
  axis.line.x = element_blank())

#==== assess sampling biases ====#
load('Data/Phanerozoic_clean_final.RData') #phanerozoic PBDB data (date accessed: 1 November 2023)
all_data <- subset(all_data, !(is.na(formation))) #remove data without formation assignments
all_data <- subset(all_data, !(formation=='')) #remove data without formation assignments

load('Data/Bioturbators_data.RData')
ecoeng_genera <- unique(bioturbators_data$genus) #get each ecosystem engineering genus name
ecoeng_formations <- unique(bioturbators_data$formation) #get each formation name containing ecosystem engineers

#Get average number of occurrences per formation for each stage
variables <- c('stage', 'n_EEforms', 'occs_EEforms', 'n_nonEEforms', 'occs_nonEEforms')
occsperform <- as.data.frame(matrix(NA, nrow=length(stage_names), ncol=length(variables)))
colnames(occsperform) <- variables

for(i in 1:nrow(occsperform)){

  this.stage <- stage_names[i]
  occsperform$stage[i] <- this.stage

  this.stage.data <- subset(all_data, stage==this.stage)

  this.stage.data_present <- subset(this.stage.data, formation %in% ecoeng_formations)
  this.stage.data_absent <- subset(this.stage.data, !(formation %in% ecoeng_formations))

  occsperform$n_EEforms[i] <- length(unique(this.stage.data_present$formation))
  occsperform$n_nonEEforms[i] <- length(unique(this.stage.data_absent$formation))

  occsperform$occs_EEforms[i] <- mean(as.data.frame(table(this.stage.data_present$formation))$Freq)
  occsperform$occs_nonEEforms[i] <- mean(as.data.frame(table(this.stage.data_absent$formation))$Freq)

}

#and same for avg number of collections per formation in each stage
variables <- c('stage', 'n_EEforms', 'colls_EEforms', 'n_nonEEforms', 'colls_nonEEforms')
collsperform <- as.data.frame(matrix(NA, nrow=length(stage_names), ncol=length(variables)))
colnames(collsperform) <- variables

for(i in 1:nrow(collsperform)){

  this.stage <- stage_names[i]
  collsperform$stage[i] <- this.stage

  this.stage.data <- subset(all_data, stage==this.stage)

  this.stage.data_present <- subset(this.stage.data, formation %in% ecoeng_formations)
  this.stage.data_absent <- subset(this.stage.data, !(formation %in% ecoeng_formations))

  collsperform$n_EEforms[i] <- length(unique(this.stage.data_present$formation))
  collsperform$n_nonEEforms[i] <- length(unique(this.stage.data_absent$formation))

  collsperform$colls_EEforms[i] <- mean(as.data.frame(table(this.stage.data_present$collection_no))$Freq)
  collsperform$colls_nonEEforms[i] <- mean(as.data.frame(table(this.stage.data_absent$collection_no))$Freq)
}

#===== load /output effect size results data  ======
load('Output/effectsizes_bioturbation_collsub.RData')
colls.form.sub_df <- bioturbation_results_df
load('Output/effectsizes_bioturbation_occsub.RData')
occs.form.sub_df <- bioturbation_results_df
load('Output/effectsizes_bioturbation_noformsub.RData')
no.form.sub_df <- bioturbation_results_df

colls.form.sub_df$sampling_difference_occs <- (occsperform$occs_EEforms - occsperform$occs_nonEEforms)
occs.form.sub_df$sampling_difference_occs <- (occsperform$occs_EEforms - occsperform$occs_nonEEforms)
no.form.sub_df$sampling_difference_occs <- (occsperform$occs_EEforms - occsperform$occs_nonEEforms)

colls.form.sub_df$sampling_difference_colls <- (collsperform$colls_EEforms - collsperform$colls_nonEEforms)
occs.form.sub_df$sampling_difference_colls <- (collsperform$colls_EEforms - collsperform$colls_nonEEforms)
no.form.sub_df$sampling_difference_colls <- (collsperform$colls_EEforms - collsperform$colls_nonEEforms)

colls.form.sub_df$method <- '5 collections per formation'
occs.form.sub_df$method  <- '20 occurrences per formation'
no.form.sub_df$method    <- 'no formation subsampling'

#order dataframes to stage sequence
colls.form.sub_df$stage <- factor(colls.form.sub_df$stage, levels=stage_names)
occs.form.sub_df$stage <- factor(occs.form.sub_df$stage, levels=stage_names)
no.form.sub_df$stage <- factor(no.form.sub_df$stage, levels=stage_names)


#===== 1 - Compare different subsampling methods ======#
#1a - Effect sizes through time 
compare_subsampling <- rbind(colls.form.sub_df, occs.form.sub_df, no.form.sub_df)
compare.cols <- c('#6aaea0', '#73628a', '#c77772')
compare.shapes <- c(21, 22, 23)

subsamp_richness <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_genrich)))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_genrich-g_genrich_sd,
                            ymax=HedgesG_genrich+g_genrich_sd,
                  fill=method),
              alpha=0.5) +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_genrich-g_genrich_sd,
                  ymax=HedgesG_genrich+g_genrich_sd,
                  color=method),
              width=5) +
  geom_line(aes(x=mid_ma, y=HedgesG_genrich, color=method)) +
  geom_point(aes(x=mid_ma, y=HedgesG_genrich, fill=method, shape=method), size=2, color='black') +
  scale_fill_manual(values=compare.cols) +
  scale_color_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_x_reverse(limits=c(538,-2), 'Time (mya)') +
  scale_y_continuous(limits=c(-1,3.2), "Hedges' g (\u00b11sd)") +
  ggtitle('Effect Size: Generic Richness') +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  sampling_theme
subsamp_richness

subsamp_H <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_H)))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                  ymax=HedgesG_H+g_H_sd,
                  fill=method),
              alpha=0.5) +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd,
                    color=method),
                width=5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, color=method)) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=method, shape=method), size=2, color='black') +
  scale_fill_manual(values=compare.cols) +
  scale_color_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_x_reverse(limits=c(538,-2), 'Time (mya)') +
  scale_y_continuous(limits=c(-1,3.2), "Hedges' g (\u00b11sd)") +
  ggtitle("Effect Size: Shannon's Diversity (H)") +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  sampling_theme
subsamp_H

subsamp_dom <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_Dominance)))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_Dominance-g_dom_sd,
                  ymax=HedgesG_Dominance+g_dom_sd,
                  fill=method),
              alpha=0.5) +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_Dominance-g_dom_sd,
                    ymax=HedgesG_Dominance+g_dom_sd,
                    color=method),
                width=5) +
  geom_line(aes(x=mid_ma, y=HedgesG_Dominance, color=method)) +
  geom_point(aes(x=mid_ma, y=HedgesG_Dominance, fill=method, shape=method), size=2, color='black') +
  scale_fill_manual(values=compare.cols) +
  scale_color_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_x_reverse(limits=c(538,-2), 'Time (mya)') +
  scale_y_continuous(limits=c(-1, 3.0), "Hedges' g (\u00b11sd)") +
  ggtitle("Effect Size: Simpson's Dominance (1/D)") +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  sampling_theme
subsamp_dom


#Effect sizes versus sampling differences
# R2 comps  - find highest between models (best fit), and then lowest between subsampling methods (worst fit)
#Richness R2 comparisons:

compare_subsampling <- subset(compare_subsampling, !is.na(sampling_difference_colls))
effort_richness_model.collsub.x1 <- lm(HedgesG_genrich~poly(sampling_difference_colls,1),
                                    data=subset(compare_subsampling, method=='5 collections per formation'))
effort_richness_model.occssub.x1 <- lm(HedgesG_genrich~poly(sampling_difference_colls,1),
                                    data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_richness_model.nosub.x1 <- lm(HedgesG_genrich~poly(sampling_difference_colls,1),
                                    data=subset(compare_subsampling, method=='no formation subsampling'))
effort_richness_model.collsub.x2 <- lm(HedgesG_genrich~poly(sampling_difference_colls,2),
                                       data=subset(compare_subsampling, method=='5 collections per formation'))
effort_richness_model.occssub.x2 <- lm(HedgesG_genrich~poly(sampling_difference_colls,2),
                                       data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_richness_model.nosub.x2 <- lm(HedgesG_genrich~poly(sampling_difference_colls,2),
                                     data=subset(compare_subsampling, method=='no formation subsampling'))
R2.richness.x1 <- c(summary(effort_richness_model.collsub.x1)$r.squared, 
                     summary(effort_richness_model.occssub.x1)$r.squared,
                     summary(effort_richness_model.nosub.x1)$r.squared)
R2.richness.x2 <- c(summary(effort_richness_model.collsub.x2)$r.squared, 
                     summary(effort_richness_model.occssub.x2)$r.squared,
                     summary(effort_richness_model.nosub.x2)$r.squared)



#Diversity R2 comparisons:
effort_diversity_model.collsub.x1 <- lm(HedgesG_H~poly(sampling_difference_colls,1),
                                       data=subset(compare_subsampling, method=='5 collections per formation'))
effort_diversity_model.occssub.x1 <- lm(HedgesG_H~poly(sampling_difference_colls,1),
                                       data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_diversity_model.nosub.x1 <- lm(HedgesG_H~poly(sampling_difference_colls,1),
                                     data=subset(compare_subsampling, method=='no formation subsampling'))
effort_diversity_model.collsub.x2 <- lm(HedgesG_H~poly(sampling_difference_colls,2),
                                       data=subset(compare_subsampling, method=='5 collections per formation'))
effort_diversity_model.occssub.x2 <- lm(HedgesG_H~poly(sampling_difference_colls,2),
                                       data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_diversity_model.nosub.x2 <- lm(HedgesG_H~poly(sampling_difference_colls,2),
                                     data=subset(compare_subsampling, method=='no formation subsampling'))
R2.diversity.x1 <- c(summary(effort_diversity_model.collsub.x1)$r.squared, 
                     summary(effort_diversity_model.occssub.x1)$r.squared,
                     summary(effort_diversity_model.nosub.x1)$r.squared)
R2.diversity.x2 <- c(summary(effort_diversity_model.collsub.x2)$r.squared, 
                     summary(effort_diversity_model.occssub.x2)$r.squared,
                     summary(effort_diversity_model.nosub.x2)$r.squared)

#Dominance R2 comparisons:
effort_dominance_model.collsub.x1 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,1),
                                        data=subset(compare_subsampling, method=='5 collections per formation'))
effort_dominance_model.occssub.x1 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,1),
                                        data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_dominance_model.nosub.x1 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,1),
                                      data=subset(compare_subsampling, method=='no formation subsampling'))
effort_dominance_model.collsub.x2 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,2),
                                        data=subset(compare_subsampling, method=='5 collections per formation'))
effort_dominance_model.occssub.x2 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,2),
                                        data=subset(compare_subsampling, method=='20 occurrences per formation'))
effort_dominance_model.nosub.x2 <- lm(HedgesG_Dominance~poly(sampling_difference_colls,2),
                                      data=subset(compare_subsampling, method=='no formation subsampling'))
R2.dominance.x1 <- c(summary(effort_dominance_model.collsub.x1)$r.squared, 
                     summary(effort_dominance_model.occssub.x1)$r.squared,
                     summary(effort_dominance_model.nosub.x1)$r.squared)
R2.dominance.x2 <- c(summary(effort_dominance_model.collsub.x2)$r.squared, 
                     summary(effort_dominance_model.occssub.x2)$r.squared,
                     summary(effort_dominance_model.nosub.x2)$r.squared)

R2_summary <- data.frame(row.names=c('5 collections/formation',
                                        '20 occurrences/formation',
                                        'no formation subsampling'),
                            R2.richness.x1, R2.richness.x2,
                            R2.diversity.x1, R2.diversity.x2,
                            R2.dominance.x1, R2.dominance.x2)
print(R2_summary)

#For y~x2:
R2_text <- as.data.frame(matrix(NA, nrow=9, ncol=3))
colnames(R2_text) <- c('measure', 'method', 'r2')
R2_text$measure <- c(rep('Richness',3),rep("Shannon's Diversity", 3),rep("Simpson's Dominance",3))
R2_text$method <- rep(unique(compare_subsampling$method),3)
R2_text$r2 <- c(sprintf("%.3f", R2_summary['5 collections/formation', 'R2.richness.x2']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.richness.x2']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.richness.x2']),
                sprintf("%.3f", R2_summary['5 collections/formation', 'R2.diversity.x2']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.diversity.x2']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.diversity.x2']),
                sprintf("%.3f", R2_summary['5 collections/formation', 'R2.dominance.x2']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.dominance.x2']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.dominance.x2'])
)
print(R2_text)

#1bii - make figures 
effort_richness <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_genrich)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_genrich,
                  fill=method, color=method), 
              method='lm', 
  #) +
              formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_genrich, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=='Richness'), x=-3, y=2.4, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Richness Effect") +
  theme_classic() +
  theme(
    strip.text=element_text(face='bold'),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    legend.position='none')

effort_diversity <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_H)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_H,
                  fill=method, color=method), 
              method='lm', 
  #) +
              formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_H, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=="Shannon's Diversity"), x=-3, y=2.1, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Shannon's Diversity Effect") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text=element_text(color='white'),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    legend.position='none')
  
effort_dominance <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_Dominance)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_Dominance,
                  fill=method, color=method), 
              method='lm', 
  #) +
              formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_Dominance, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=="Simpson's Dominance"), x=-3, y=2.7, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Simpson's Dominance Effect") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text=element_text(color='white'),
    legend.title=element_blank(),
    legend.position='none')

#Compare all together
effort_comps <- ggarrange2(effort_richness, effort_diversity, 
          ncol=1)


#for y~x:
R2_text <- as.data.frame(matrix(NA, nrow=9, ncol=3))
colnames(R2_text) <- c('measure', 'method', 'r2')
R2_text$measure <- c(rep('Richness',3),rep("Shannon's Diversity", 3),rep("Simpson's Dominance",3))
R2_text$method <- rep(unique(compare_subsampling$method),3)
R2_text$r2 <- c(sprintf("%.3f", R2_summary['5 collections/formation', 'R2.richness.x1']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.richness.x1']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.richness.x1']),
                sprintf("%.3f", R2_summary['5 collections/formation', 'R2.diversity.x1']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.diversity.x1']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.diversity.x1']),
                sprintf("%.3f", R2_summary['5 collections/formation', 'R2.dominance.x1']),
                sprintf("%.3f", R2_summary['20 occurrences/formation', 'R2.dominance.x1']),
                sprintf("%.3f", R2_summary['no formation subsampling', 'R2.dominance.x1'])
)
print(R2_text)

#1bii - make figures 
effort_richness <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_genrich)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_genrich,
                  fill=method, color=method), 
              method='lm', 
  ) +
  #formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_genrich, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=='Richness'), x=-3, y=2.4, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Richness Effect") +
  theme_classic() +
  theme(
    strip.text=element_text(face='bold'),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    legend.position='none')

effort_diversity <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_H)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_H,
                  fill=method, color=method), 
              method='lm', 
  ) +
  #formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_H, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=="Shannon's Diversity"), x=-3, y=2.1, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Shannon's Diversity Effect") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text=element_text(color='white'),
    legend.title=element_blank(),
    axis.title.x=element_blank(),
    legend.position='none')

effort_dominance <- ggplot(data=subset(compare_subsampling, !(is.na(HedgesG_Dominance)))) +
  geom_smooth(aes(x=sampling_difference_colls, y=HedgesG_Dominance,
                  fill=method, color=method), 
              method='lm', 
  ) +
  #formula=y~poly(x,2)) +
  geom_point(aes(x=sampling_difference_colls, y=HedgesG_Dominance, 
                 fill=method, shape=method), size=2, colour='black') +
  geom_text(data=subset(R2_text, measure=="Simpson's Dominance"), x=-3, y=2.7, hjust=0, aes(label=paste0("R^2~'='~", r2)), parse = TRUE) +
  facet_wrap(vars(method)) +
  scale_fill_manual(values=compare.cols) +
  scale_shape_manual(values=compare.shapes) +
  scale_color_manual(values=compare.cols) +
  scale_y_continuous("Hedges' g") +
  scale_x_continuous("avg. n. collections per EE formations — avg. n. collections per non EE formations") +
  ggtitle("Sampling Effort versus Simpson's Dominance Effect") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text=element_text(color='white'),
    legend.title=element_blank(),
    legend.position='none')

effort_comps <- ggarrange2(effort_richness, effort_diversity, 
                           ncol=1)
