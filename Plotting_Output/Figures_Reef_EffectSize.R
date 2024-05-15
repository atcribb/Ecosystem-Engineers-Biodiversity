#Authors: Alison Cribb, Will Gearty
#Summary: Final analyses and output plots for reef effect sizes 

# clear old data
rm(list = ls())

#==== load packages ====#
library(divDyn)
data("stages", package="divDyn")
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]
period_names <- unique(stages[which(stages$stage %in% stage_names), 'system'])
period.cols <- unique(stages[which(stages$stage %in% stage_names), 'systemCol'])
library(ggplot2)
library(deeptime)
library(metafor)
library(tidyr)
library(RColorBrewer)

effectsize_theme <- theme(
  legend.position="inside",
  legend.position.inside=c(0.9, 0.9),
  legend.title=element_blank(),
  plot.title=element_text(hjust=0.5),
  axis.text = element_text(color = "black"),
  axis.line.x = element_blank())

#==== load data ====#
load('Output/effectsizes_reefs_occsub.RData')
results_df <- reef_results_df

results_df$stage <- factor(results_df$stage, levels=stage_names)
results_df$period <- factor(results_df$period, levels=period_names)

#===========  MAIN EFFECT SIZES FIGURES ===========

#Richness effect size
#restructure M1 vs M2 data
m1_richness <- results_df[,c(1:3,6:7)]
colnames(m1_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m1_richness$formations <- 'Reef-builders present'
m2_richness <- results_df[,c(1:3,8:9)]
colnames(m2_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m2_richness$formations <- 'Reef-builders absent'
compare_richness <- rbind(m1_richness, m2_richness)

# these are pretty similar for colorblind folks, might want to change one of them -WG
compare.forms.cols <- c('#A03544', '#B19398')
compare_richness$formations <- factor(compare_richness$formations, 
                                      levels=c('Reef-builders present', 'Reef-builders absent'))

compare_richness_plot <- ggplot(data=compare_richness) +
  #geom_ribbon(aes(x=mid_ma, ymin=richness-sd, ymax=richness+sd, fill=formations), alpha=0.5) +
  geom_errorbar(aes(x=mid_ma, ymin=richness-sd, ymax=richness+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=richness, colour=formations)) +
  geom_point(aes(x=mid_ma, y=richness, fill=formations), colour='black', size=3, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,15), name="Generic Richness") +
  ggtitle('Comparison of Generic Richness') +
  coord_geo(pos="bottom", dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
compare_richness_plot  

results_df$genrich_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$HedgesG_genrich) < 0.2),'genrich_significance'] <- 'no effect'
results_df[which( abs(results_df$HedgesG_genrich) >= 0.2 & abs(results_df$HedgesG_genrich) < 0.5),'genrich_significance'] <- 'weak effect'
results_df[which( abs(results_df$HedgesG_genrich) >= 0.5 & abs(results_df$HedgesG_genrich) < 0.8),'genrich_significance'] <- 'medium effect'
results_df[which( abs(results_df$HedgesG_genrich) >= 0.8),'genrich_significance'] <- 'strong effect'
results_df$genrich_significance <- factor(results_df$genrich_significance, levels=c('no effect', 'weak effect', 'medium effect', 'strong effect'))


effectsize_richness_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_genrich))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_genrich-g_genrich_sd,
                    ymax=HedgesG_genrich+g_genrich_sd),
                width=8) +
  geom_point(aes(x=mid_ma, y=HedgesG_genrich, size=genrich_significance), fill='white', shape=21) +
  geom_point(aes(x=mid_ma, y=HedgesG_genrich, fill=period, alpha=genrich_significance, size=genrich_significance), shape=21) +
  scale_size_manual(values=c(2, 2.4, 3.2, 3.6)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_fill_manual(values=period.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,2.5), name="Hedges' g (\u00b11sd)") +
  guides(fill='none') +
  ggtitle('Reef-builder Effect Size: Generic Richness') +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
effectsize_richness_plot  

full_richness_fig <- ggarrange2(compare_richness_plot, effectsize_richness_plot, ncol=1)

#Shannon's Diversity
#restructure M1 vs M2 data
m1_H <- results_df[,c(1:3,12:13)]
colnames(m1_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m1_H$formations <- 'Reef-builders present'
m2_H <- results_df[,c(1:3,14:15)]
colnames(m2_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m2_H$formations <- 'Reef-builders absent'
compare_H <- rbind(m1_H, m2_H)

compare.forms.cols <- c('#A03544', '#B19398')
compare_H$formations <- factor(compare_H$formations, 
                               levels=c('Reef-builders present', 'Reef-builders absent'))

compare_H_plot <- ggplot(data=compare_H) +
  geom_errorbar(aes(x=mid_ma, ymin=H-sd, ymax=H+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=H, colour=formations)) +
  geom_point(aes(x=mid_ma, y=H, fill=formations), colour='black', size=3, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,2.5), name="Shannon's Diversity (H)") +
  ggtitle("Comparison of Shannon's Diversity") +
  coord_geo(pos="bottom", dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
compare_H_plot  

results_df$H_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
results_df[which( abs(results_df$HedgesG_H) >= 0.2 & abs(results_df$HedgesG_H) < 0.5),'H_significance'] <- 'weak effect'
results_df[which( abs(results_df$HedgesG_H) >= 0.5 & abs(results_df$HedgesG_H) < 0.8),'H_significance'] <- 'medium effect'
results_df[which( abs(results_df$HedgesG_H) >= 0.8),'H_significance'] <- 'strong effect'
results_df$H_significance <- factor(results_df$H_significance, levels=c('no effect', 'weak effect', 'medium effect' , 'strong effect'))

effectsize_H_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_H))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=8) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance), fill='white', shape=21) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=period, alpha=H_significance, size=H_significance), shape=21) +
  scale_size_manual(values=c(2, 2.4, 3.2, 3.6)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_fill_manual(values=period.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1,2.3), name="Hedges' g (\u00b11sd)") +
  guides(fill='none') +
  ggtitle('Reef-builder Effect Size: Shannons Diversity') +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
effectsize_H_plot  

full_H_fig <- ggarrange2(compare_H_plot, effectsize_H_plot, ncol=1)

#Simpson's Dominance
#restructure M1 vs M2 data
m1_dom <- results_df[,c(1:3,18:19)]
colnames(m1_dom) <- c('period', 'stage', 'mid_ma', 'dom', 'sd')
m1_dom$formations <- 'Reef-builders present'
m2_dom <- results_df[,c(1:3,20:21)]
colnames(m2_dom) <- c('period', 'stage', 'mid_ma', 'dom', 'sd')
m2_dom$formations <- 'Reef-builders absent'
compare_dom <- rbind(m1_dom, m2_dom)

compare.forms.cols <- c('#A03544', '#B19398')
compare_dom$formations <- factor(compare_dom$formations, 
                                 levels=c('Reef-builders present', 'Reef-builders absent'))

compare_dom_plot <- ggplot(data=compare_dom) +
  geom_errorbar(aes(x=mid_ma, ymin=dom-sd, ymax=dom+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=dom, colour=formations)) +
  geom_point(aes(x=mid_ma, y=dom, fill=formations), colour='black', size=3, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(name="Simpson's Dominance (1/D)") +
  ggtitle("Comparison of Simpson's Dominance") +
  coord_geo(pos="bottom", dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
compare_dom_plot 

results_df$dominance_significance <- rep(NA, nrow(results_df))
# 0.2	=>	small effect
# 0.5	=>	medium effect
# 0.8	=>	large effect

results_df$dominance_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$HedgesG_Dominance) < 0.2),'dominance_significance'] <- 'no effect'
results_df[which(abs(results_df$HedgesG_Dominance) >= 0.2 & abs(results_df$HedgesG_Dominance) < 0.5),'dominance_significance'] <- 'weak effect'
results_df[which(abs(results_df$HedgesG_Dominance) >= 0.5 & abs(results_df$HedgesG_Dominance) < 0.8),'dominance_significance'] <- 'medium effect'
results_df[which(abs(results_df$HedgesG_Dominance) >= 0.8),'dominance_significance'] <- 'strong effect'
results_df$dominance_significance <- factor(results_df$dominance_significance, levels=c('no effect', 'weak effect', 'medium effect', 'strong effect'))

effectsize_dom_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_Dominance))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_Dominance-g_dom_sd,
                    ymax=HedgesG_Dominance+g_dom_sd),
                width=8) +
  geom_point(aes(x=mid_ma, y=HedgesG_Dominance, size=dominance_significance), fill='white', shape=21) +
  geom_point(aes(x=mid_ma, y=HedgesG_Dominance, fill=period, alpha=dominance_significance, size=dominance_significance), shape=21) +
  scale_fill_manual(values=period.cols) +
  scale_size_manual(values=c(2, 2.4, 3.2, 3.6)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1,2.5), name="Hedges' g (\u00b11sd)") +
  guides(fill='none') +
  ggtitle("Reef-builder Effect Size: Simpson's Dominance") +
  coord_geo(pos='bottom', dat='periods', size='auto', abbrv=FALSE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
effectsize_dom_plot  

full_dominance_fig <- ggarrange2(compare_dom_plot, effectsize_dom_plot, ncol=1)
