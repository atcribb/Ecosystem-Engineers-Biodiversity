# Marine Phanerozoic biodiversity increases in presence of ecosystem engineers

<i><b> THIS README IS UNDER CONSTRUCTION </i></b>

If you download this repository, set it as your working directory and everything should run smoothly. You will need to create an /Output folder to save the results .RData files produced in the analyses so that you can then generate the figures.

## Datasets
There are three .RData files in <kbd>[Data](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Data)</kbd>:
* ``Bioturbators_data.RData`` - contains the bioturbator ecosystem engineering data
* ``Reef_Ecosystem_Engineers_Final.RData`` - contains the reef-builder ecosystem engineering data
* ``Phanerozoic_clean_final.Rdata`` - contains the cleaned and stratigraphically binned (largely following Kocsis et al. 2019 ddPhanero protocol with additional cleaning of formation data) Phanerozoic data. The raw PBDB data was downloaded 1 November 2023. Using this dataset will reproduce the figures in the paper. 

Analyases are not done on these two datasets perse, but they are needed in order to come up with the lists of formations that do and do not contain ecosystem engineers in the <b>Analyses</b> Effect Size scripts.


## Analyses
There are two types of <kbd>[Analyses](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Data)</kbd>: the main effect size analyses, and the sampling biases analyses.

* ``EffectSize_Bioturbation.R`` 
* ``EffectSize_Reefs.R``


## Plotting outputs 

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

