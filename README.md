# Cross_Environment_Global_Epistasis
This repository contains all relevant code for completing analyses and producing figures for the paper 'Environment-independent distribution of mutational effects emerges from microscopic epistasis'. The zenodo DOI associated is 10.5281/zenodo.12696536. To use this code, you will need an updated version of R.  To properly run this code ensure:
1. You have installed all packages listed at the top of the script.
2. You update the working directory (and, optionally, the path to save files) in the initial lines of the code. The code will only run if your working directory has all relevant files.

Unfortunately, the confirmation experiment BC count file is too large to upload here. It is available upon request. Please note that much of the script associated with the Reddy-Desai model (E.g., supplementary figures S12-13) are directly reproduced with permission from Reddy and Desai eLife 2021 (https://github.com/greddy992/global_epistasis). 

The file “CoreFigureCode.R”, when run with the relevant data files in the working directory, will produce all main text and supplemental figures. Note that the following files are required and readily available in this repository):
-	**df_s_ests_wGR** : Contains the average fitness effect estimates for all mutations in all strains and environments where we had sufficient measurements. Column names are as follows
-	  Strain: Strain names, consistent with Johnson et al 2019
-	  Env: environment name, as ‘Temp’ ‘SC’ ‘pH’
-	  Mut_ID: Numerical identifier corresponding to each mutation. Gene names for each can be found in Mut_fullinfo_use data frame
-	  avgS: Average fitness effect of each mutation
-	  varS: variance in the fitness effect
-	  seS: standard error of the fitness effect
-	  Mean_GR: Estimated growth rate for the corresponding strain in the corresponding environment
-	  Std_err: Standard error of the growth rate estimate
-	  mutCall: call of each mutation’s effect
-	  minCI: bottom of the 99% confidence interval range for the mutation’s effect
-	  maxCI: top of the 99% confidence interval range for the mutation’s effect
-	**growthRate_data** : Contains the estimated growth rate all all background strains in all environments
-	**elife-27167-supp1-v2.csv**: Obtained from Jerison et al eLife 2017. Contains indicators of the allelic variant (0 or 1; BY vs RM) at each focal QTL locus (indicated as column names)
-	**Bloom_QTLS_tempPH.csv** : A subset of data from Bloom et al Nature 2013. Contains all temperature and pH relevant QTLS annotated in their study.
-	**Mut_fullInfo_use**: Maps mutation numerical IDs to the gene they are in or nearby
-	**JohnsonDFEmeanDat** – can be created given data from Johnson et al Science 2019: simply need the strain names and the corresponding parameters: DFEmean, DFEvar, DFEskew , DFEmean_min, DFEmean_max, DFEvar_min, DFEvar_max, DFEskew_min, DFEskew_max. I have pre-processed this freely available data and uploaded the file here for ease of use.
-	  However: The growth rates provided in Johnson et al 2019 are relative growth rate metrics. We convert these to exponential growth rate using the relationship between relative and exponential growth given in their supplemental figure S9 and detailed in the supplementary materials of the manuscript accompanying this code.
-	**allExpGR_data.csv**: Contains estimates of the raw exponential growth rate of all mutants in the study
-	**Milo_S_est.csv**: Contains fitness effects of mutations estimated from Johnson et al Science 2019
-	**milo_slopes_ints**: A processed version of the above file with the best fit global epistasis slopes and intercepts from Johnson et al 2019
-	**allBC_s_ests.csv**: This contains the estimated fitness effects for all barcodes, and is neccessary only for Figure S3C. Unfortunately, it is to large to uplaod here, but it is available on request - please reach out to the corresponding author of the manuscript


Note that many sections depend on data files created in other sections, so it is best run all together, Total runtime is ~4 hours (on Mac 16GB Ram M1 chip laptop), but the majority of this is in the simulations in the last section (Figure S13). Total runtime in ~15 minutes without S13. 
