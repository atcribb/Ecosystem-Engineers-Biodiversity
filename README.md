# Marine Phanerozoic biodiversity increases in presence of ecosystem engineers

<i><b> THIS README IS UNDER CONSTRUCTION </i></b>

#Datasets
There are two .RData files in <b>Data</b>:
*Bioturbators_data.RData - contains the bioturbator ecosystem engineering data
*Reef_Ecosystem_Engineers_Final.RData - contains the reef-builder ecosystem engineering data

Analyases are not done on these two datasets perse, but they are needed in order to come up with the lists of formations that do and do not contain ecosystem engineers in the <b>Analyses</b> Effect Size scripts.

#Data Cleaning 
These scripts will help you assemble a cleaned and stage-binned PBDB database. The following steps will create a <b>Phanerozoic_clean_final.RData</b> dataset that is used in the analyses. Unfortunately, the PBDB cannot be uploaded to github... 
*Step 1: Download marine PBDB data. You can use the PBDB URL in <b>pbdb_cleaning_functions.R</b>
*Step 2: Use the clean_marine function in <b>pbdb_cleaning_functions</b> on your raw PBDB download. Save this with some appropriate name (ie pbdb_cleaned_notbinned.RData). This cleaning function follows https://github.com/divDyn/ddPhanero, Kocsis et al. 2018 with the divDyn package
*Step 3: Use the stage_binning function in <b>pbdb_binning_functions</b> on your cleaned PBDB data. Save this as Phanerozoic_marine_cleaned_binned.RData. This binning function follows https://github.com/divDyn/ddPhanero, Kocsis et al. 2018 with the divDyn package. 
*Step 4: The analyses are conducted by comparing ecological metrics between formations, but your PBDB data at this stage has an issue with formation synonyms, typos, language translations, incorrect information, and more. This means that fossils that actually exist together in one formation are entered into the PBDB in many different formations, and we need to fix this before we can run our analyses. Use the <b>formation_cleaning.R</b> script and the <b>formation_sorting.csv</b> file to clean up the PBDB data at the formation level so that the data appropriately falls into ecosystem engineering versus non-ecosystem engineering bins. Save this final data as Phanerozoic_clean_final.RData


#Analyses

#Plotting outputs 
