# Marine Phanerozoic biodiversity increases in presence of ecosystem engineers

<i><b> THIS README IS UNDER CONSTRUCTION </i></b>

## Datasets
There are three .RData files in <kbd>[Data](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Data)</kbd>:
* ``Bioturbators_data.RData`` - contains the bioturbator ecosystem engineering data
* ``Reef_Ecosystem_Engineers_Final.RData`` - contains the reef-builder ecosystem engineering data
* ``Phanerozoic_clean_final.Rdata`` - contains the cleaned and stratigraphically binned (largely following Kocsis et al. 2019 ddPhanero protocol with additional cleaning of formation data) Phanerozoic data. The raw PBDB data was downloaded 1 November 2023. Using this dataset will reproduce the figures in the paper. 

 ``Bioturbators_data.RData`` and ``Reef_Ecosystem_Engineers_Final.RData`` are used to identify ecosystem engineering genera and formations that contain ecosystem engineering genera. Actual effect sizes analyses are performed on Phanerozoic_clean_final.RData 

## Analyses
There are two sets of analyses in <kbd>[Analyses](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Analyses)</kbd>: The main effect size analyses, and sampling biases analyses 

### Effect Size Analyses: 
<b> Please be aware that these scripts can take more than 24 hours to run with the 1000 iterations used in the paper's final results. Consider decreasing number of iterations if using personal computer. </b>

* ``EffectSize_Bioturbators.R`` - calculate Hedges g effect size values for bioturbators impacts on biodiversity
* ``EffectSize_Bioturbators.R`` - calculates Hedges g effect size values for reef-builder impacts on biodiversity
* <b> These analyses need to be run first to create output files. </b> You will need to create an /Output folder in your working directory to save the results. There are three subsampling protocols to choose from at the top of the script that change how formations are subsampled. The subsampling protocol corresponds to file save line at the end of the script. For example, if subsampling 20 occurrences per formation (as used in the final results), leave form.subsampling <- 'occurrences' and save(bioturbation_results_df, file='Output/effectsizes_bioturbation_occsub.RData') not commented out. You do not need to change anything else in the script other than the subsampling selection and appropriate file save. 

### Sampling Biases Analyses:
* ``Bioturbation_SamplingBiases.R`` - determines which formation subsampling protocol best minimizes sampling biases on bioturbation effect sizes.
* ``Reefs_SamplingBiases.R`` - determines which formations subsampling protocol best minimizes sampling biases on reef-builder effect sizes 
* These analyses will call results from all three subsampling protocols in the effect size analyses. To run the sampling biases analyses, be sure all results are saved to an /Output file.

## Plotting outputs 
* ``Figures_Bioturbation_EffectSizes.R`` - creates the main text figures for the effects of bioturbation on biodiversity indices.
* ``Figures_Reef_EffectSize.R`` - creates the main text figures for the effects of reef-buidlers on biodiversity indices.
* These by default load in the results from subsampling 20 occurrences per formation to reproduce the figures in the paper. You may change these to any of your results files in /Output generated from the Effect Size Analyses.

# Change log
### 25 April 2024
```diff
+ upload full PBBD dataset
+ update README to reflect new dataset upload
+ update README for Analyses and Plotting_Output
- remove DataClean in lieu of uploading PBDB data
```

### 4 April 2024
```diff
+ writing the README
+ uploading Analyses scripts
```

