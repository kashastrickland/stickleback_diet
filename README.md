This repo contains all data and code required to reproduce results, tables and figures presented in the manuscript entitled "Spatial variation and individual specialization of stickleback diet 
in relation to trophic morphology" by Snorradóttir et al..

Data are contained within the /data subfolder. Contained in the folder are the following main dataframes:

"**dissection_data.csv**" contains the raw format of diet and morphological data used in analyses  
"**Myvatn_WSGUTM28.shp**" is the shape file require to reproduce the map in figure 1  
Files ending in "**.TPS**" are the landmarked data required for geometric morphometric analyses

R code to reproduce all results include:

"**Head_morphology.R** which is required to run geometric morphometric analyses and reproduce figure 3. 
"**Map Figure.R** reproduces Figure 1A
"**Prepare Data.R** organises the raw data for analysis as described in the manuscript, and also reproduced Figure 1B and table S1.
"**morphdivergence_And_dietspecialisation.R** reproduced RDA analysis, analysis for spatial divergence in morphology, analyses of PSi (i.e., diet specialisation). Also reproduces tables 1, 3 and 4, and figures 3 and 6.


