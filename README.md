This repo contains all data and code required to reproduce results, tables and figures presented in the manuscript entitled "_Spatial variation and individual specialization of stickleback diet 
in relation to trophic morphology_" by Snorradóttir et al..

Data are contained within the **/data** subfolder. Contained in the folder are the following main dataframes:

"**dissection_data.csv**" contains the raw format of diet and morphological data used in analyses  
"**Myvatn_WSGUTM28.shp**" is the shape file require to reproduce the map in figure 1  
Files ending in "**.TPS**" are the landmarked data required for geometric morphometric analyses  

R code to reproduce all results include:

"**Head_morphology.R**" which is required to run geometric morphometric analyses and reproduce figure 3.  
"**Map Figure.R**" reproduces Figure 1A.  
"**Prepare Data.R**" organises the raw data for analysis as described in the manuscript, and also reproduced Figure 1B and table S1.  
"**morphdivergence_And_dietspecialisation.R**" reproduced RDA analysis, analysis for spatial divergence in morphology, analyses of PSi (i.e., diet specialisation). Also reproduces tables 1, 3 and 4, and figures 3 and 6.   
"**Fit Taxon-Sppecific Models.R**" runs models for diet variation through space and according to morphology.  
"**Derived Quantities for Visualisation.R**" processes model output from Fit Taxon-Sppecific Models.R to prepare for visualisation in Figures 4 and 5.  
"**Visualise Taxon-Specific Responses.R**" reproces Figures 4 and 5.  

  For all code explained above authorship is noted at the top of the code.  


  All figures and tables in **/tables_figures**, including inkscape (graphical software) files used to produce multipanel plots.



