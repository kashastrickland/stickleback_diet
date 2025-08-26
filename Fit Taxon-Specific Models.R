#=========================================================================================
#========== Preliminaries
#=========================================================================================

# Load packages
library(tidyverse)
library(brms)

# Set up parallel processing (for brms)
options(mc.cores = parallel::detectCores()-2)

# import data
diet_long_analysis <- read_csv("data/diet_long_analysis.csv")

#=========================================================================================



#=========================================================================================
#========== Fit models
#=========================================================================================

# fit full model
# m_diet_full <- brm(bf(count_in_stomach ~ length_z + gut_length_z + Gap_z + raker_length_z +
#                         Comp1_z + Comp2_z + Comp3_z + sex + 
#                         (1 | taxon * Site + taxon:sex + fish_id) +
#                         (0 + length_z | taxon) +
#                         (0 + gut_length_z | taxon) +
#                         (0 + Gap_z | taxon) +
#                         (0 + raker_length_z | taxon) +
#                         (0 + Comp1_z | taxon) +
#                         (0 + Comp2_z | taxon) +
#                         (0 + Comp3_z | taxon),
#                       zi ~ (1 | taxon)),
#                    data = diet_long_analysis,
#                    family = zero_inflated_poisson("log"),
#                    iter = 4000,
#                    chains = 4,
#                    cores = 4,
#                    control = list(adapt_delta = 0.99))
# write_rds(m_diet_full, "model_fits/m_diet_full.rds")

#=========================================================================================





