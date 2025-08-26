#############################
## Author: Joseph Phillips, 2025

#=========================================================================================
#========== Preliminaries
#=========================================================================================

# Load packages
library(tidyverse)
library(brms)

# Set up parallel processing 
options(mc.cores = parallel::detectCores()-2)

# import data
diet_long_analysis <- read_csv("data/diet_long_analysis.csv")

# import model fit
m_diet <- read_rds("model_fits/m_diet_full.rds")

# define taxon levels and labels
taxon_levels <- c("adult_fly", 
                  "alona", 
                  "chrinominae", 
                  "copepods", 
                  "mollusks", 
                  "orthocladinae",     
                  "ostracod",          
                  "other_cladocerans", 
                  "stickleback_eggs",  
                  "tanypodinae")
taxon_labels <- c("Chironomidae (adult)", 
                  "Alona spp.", 
                  "Chrinominae", 
                  "Cyclopideae", 
                  "Mollusca", 
                  "Orthocladinae",     
                  "Ostracoda",          
                  "Cladocerans (other)", 
                  "Stickleback eggs",  
                  "Tanypodinae")

#=========================================================================================





#=========================================================================================
#========== SDs characterizing taxon-specific responses
#=========================================================================================

# extract parameters and convert to long form
# assign tidy names
slope_sd_taxon <- tibble(as.data.frame(m_diet, 
                                  variable = c("sd_taxon__length_z",
                                               "sd_taxon__gut_length_z",
                                               "sd_taxon__Gap_z",
                                               "sd_taxon__raker_length_z",
                                               "sd_taxon__Comp1_z",
                                               "sd_taxon__Comp2_z",
                                               "sd_taxon__Comp3_z"))) %>%
  pivot_longer(1:7) %>%
  mutate(name = factor(name, 
                       levels = c("sd_taxon__length_z",
                                  "sd_taxon__gut_length_z",
                                  "sd_taxon__Gap_z",
                                  "sd_taxon__raker_length_z",
                                  "sd_taxon__Comp1_z",
                                  "sd_taxon__Comp2_z",
                                  "sd_taxon__Comp3_z"),
                       labels = c("Body Length",
                                  "Gut Length",
                                  "Gap",
                                  "Raker Length",
                                  "PC1",
                                  "PC2",
                                  "PC3")))

# export
# write_csv(slope_sd_taxon, "derived_quantities/slope_sd_taxon.csv")

#=========================================================================================




#=========================================================================================
#========== Taxon-specific slopes
#=========================================================================================

# define object with taxon-specific coefficients 
# summarize using intervals 68% (SEs) and 95% coverage
taxon_coefs <- coef(m_diet, probs = c(0.05, 0.16, 0.84, 0.95))$taxon 

# extract coefficient names
coef_names <- c("length_z", "gut_length_z", "Gap_z", "raker_length_z",
                "Comp1_z", "Comp2_z", "Comp3_z")

# reformat into long data frame
# assign tidy taxon names
taxon_coefs_long <- lapply(coef_names, function(x_){
  y_ = taxon_coefs[, , x_]
  tibble(as.data.frame(y_)) %>% mutate(coef = x_, taxon = rownames(y_))
}) %>% 
  bind_rows() %>%
  select(taxon, coef, Estimate, Est.Error, Q5, Q16, Q84, Q95) %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# export
# write_csv(taxon_coefs_long, "derived_quantities/taxon_coefs_long.csv")

#=========================================================================================




#=========================================================================================
#========== SDs characterizing intercepts
#=========================================================================================

# extract parameters and convert to long form
# assign tidy names
intercept_sd <- tibble(as.data.frame(m_diet, 
                                       variable = c("sd_fish_id__Intercept",
                                                    "sd_Site__Intercept",
                                                    "sd_taxon__Intercept",
                                                    "sd_taxon:Site__Intercept",
                                                    "sd_taxon:sex__Intercept"))) %>%
  pivot_longer(1:5) %>%
  mutate(name = factor(name, 
                       levels = c("sd_fish_id__Intercept",
                                  "sd_Site__Intercept",
                                  "sd_taxon__Intercept",
                                  "sd_taxon:Site__Intercept",
                                  "sd_taxon:sex__Intercept"),
                       labels = c("Fish ID",
                                  "Site",
                                  "Taxon",
                                  "Taxon x Site",
                                  "Taxon x Sex")))

# export
# write_csv(intercept_sd, "derived_quantities/intercept_sd.csv")

#=========================================================================================





#=========================================================================================
#========== Spatial variation in taxonomic representation in diet
#=========================================================================================

# extract taxon variation among sites
taxon_intercepts_sites_raw <- tibble(rownames_to_column(
  as.data.frame(coef(m_diet, probs = c(0.05, 0.16, 0.84, 0.95))$`taxon:Site`[ , , 1])))

# define function to extract taxon names
taxon_f <- function(x_){
  x1_ = strsplit(x_, "_")[[1]]
  i = length(x1_)
  x2_ <- x1_[1:(i - 1)]
  x3_ <- paste(x2_, collapse = "_")
  return(x3_)
}

# define function to extract site names
site_f <- function(x_){
  x1_ = strsplit(x_, "_")[[1]]
  i = length(x1_)
  x2_ <- x1_[i]
  return(x2_)
}

# create clean data with taxon and site names
taxon_intercepts_sites <- taxon_intercepts_sites_raw %>% 
  mutate(taxon = rowname %>% map_chr(~taxon_f(.x)),
         Site = rowname %>% map_chr(~site_f(.x))) %>%
  select(taxon, Site, Estimate, Est.Error, Q5, Q95, Q16, Q84) %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# export
# write_csv(taxon_intercepts_sites, "derived_quantities/taxon_intercepts_sites.csv")

#=========================================================================================





#=========================================================================================
#========== Spatial variation among taxon intercepts
#=========================================================================================

# extract taxon intercepts
taxon_int <- tibble(rownames_to_column(
  as.data.frame(coef(m_diet, probs = c(0.05, 0.16, 0.84, 0.95))$`taxon`[ , , 1])))%>%
  rename(taxon = rowname) %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))


# export
# write_csv(taxon_int, "derived_quantities/taxon_int.csv")

#=========================================================================================







#=========================================================================================
#========== Sex differences in taxonomic representation in diet
#=========================================================================================

taxon_sex_comb <- as.vector(outer(taxon_levels, 
                                  unique(diet_long_analysis$sex), 
                                  paste, 
                                  sep="_"))
sex_ranef_vars <- tibble(taxon_c = paste("r_taxon:sex[", 
                                         taxon_sex_comb, 
                                         ",Intercept]", 
                                         sep=""),
                         taxon = rep(taxon_levels, 2),
                         sex = rep(unique(diet_long_analysis$sex), each = length(taxon_levels)))

taxon_interepts_F_relto_M <- as.data.frame(m_diet,
              variable = sex_ranef_vars$taxon_c) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(1:length(sex_ranef_vars$taxon_c),
               names_to = "taxon_c") %>%
  # add tidy taxon labels
  full_join(sex_ranef_vars) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  pivot_wider(names_from = sex, values_from = value) %>%
  mutate(value = `F` - M) %>%
  group_by(taxon) %>%
  summarize(Estimate = quantile(value, probs = 0.5),
            Q5 = quantile(value, probs = 0.025),
            Q16 = quantile(value, probs = 0.16),
            Q84 = quantile(value, probs = 0.84),
            Q95 = quantile(value, probs = 0.975)) %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))
  
# write_csv(taxon_interepts_F_relto_M, "derived_quantities/taxon_interepts_F_relto_M.csv")

#=========================================================================================





#=========================================================================================
#========== Zero-Inflation by Taxon
#=========================================================================================

# extract zero inflated values
zis <- ranef(m_diet, summary = FALSE)$taxon[,,"zi_Intercept"] + 
  fixef(m_diet, summary = FALSE)[, "zi_Intercept"]

# convert linear predictor to zi probability
# convert to long scale
# summarize by taxon
taxon_zi <- pivot_longer(as.data.frame(exp(zis) / (1 + exp(zis))), 1:10, names_to = "taxon") %>%
  group_by(taxon) %>%
  summarize(Estimate = quantile(value, probs = 0.5),
            Q5 = quantile(value, probs = 0.025),
            Q16 = quantile(value, probs = 0.16),
            Q84 = quantile(value, probs = 0.84),
            Q95 = quantile(value, probs = 0.975)) %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(taxon_zi, "derived_quantities/taxon_zi.csv")

#=========================================================================================





#=========================================================================================
#========== Fitted values: Length
#=========================================================================================

fitted_length <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_length_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",length_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                                       ",length_z]", 
                                       sep=""),
                       taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_length_z,
                 r_b,
                 r_zi),
         length_z = seq(min(diet_long_analysis$length_z),
                        max(diet_long_analysis$length_z),
                        length.out = 10)) %>%
  # calculate linear predictor
  # inverse logit transform to calculate zi
  # calculate expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_length_z + r_b) * length_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           length_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # convert predictor variable to natural scale
  mutate(length = 
           sd(diet_long_analysis$length) * length_z + mean(diet_long_analysis$length)) %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_length, "derived_quantities/fitted_length.csv")


#=========================================================================================





#=========================================================================================
#========== Fitted values: Gut Length
#=========================================================================================

fitted_gut_length <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_gut_length_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",gut_length_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",gut_length_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_gut_length_z,
                 r_b,
                 r_zi),
         gut_length_z = seq(min(diet_long_analysis$gut_length_z),
                        max(diet_long_analysis$gut_length_z),
                        length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_gut_length_z + r_b) * gut_length_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           gut_length_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_gut_length, "derived_quantities/fitted_gut_length.csv")


#=========================================================================================





#=========================================================================================
#========== Fitted values: Gap
#=========================================================================================

fitted_gap <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_Gap_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",Gap_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",Gap_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_Gap_z,
                 r_b,
                 r_zi),
         Gap_z = seq(min(diet_long_analysis$Gap_z),
                            max(diet_long_analysis$Gap_z),
                            length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_Gap_z + r_b) * Gap_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           Gap_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_gap, "derived_quantities/fitted_gap.csv")


#=========================================================================================





#=========================================================================================
#========== Fitted values: Raker Length
#=========================================================================================

fitted_raker_length <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_raker_length_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",raker_length_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",raker_length_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_raker_length_z,
                 r_b,
                 r_zi),
         raker_length_z = seq(min(diet_long_analysis$raker_length_z),
                            max(diet_long_analysis$raker_length_z),
                            length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_raker_length_z + r_b) * raker_length_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           raker_length_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_raker_length, "derived_quantities/fitted_raker_length.csv")

#=========================================================================================





#=========================================================================================
#========== Fitted values: Comp1
#=========================================================================================

fitted_Comp1 <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_Comp1_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",Comp1_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",Comp1_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_Comp1_z,
                 r_b,
                 r_zi),
         Comp1_z = seq(min(diet_long_analysis$Comp1_z),
                              max(diet_long_analysis$Comp1_z),
                              length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_Comp1_z + r_b) * Comp1_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           Comp1_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_Comp1, "derived_quantities/fitted_Comp1.csv")

#=========================================================================================





#=========================================================================================
#========== Fitted values: Comp2
#=========================================================================================

fitted_Comp2 <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_Comp2_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",Comp2_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",Comp2_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_Comp2_z,
                 r_b,
                 r_zi),
         Comp2_z = seq(min(diet_long_analysis$Comp2_z),
                              max(diet_long_analysis$Comp2_z),
                              length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_Comp2_z + r_b) * Comp2_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           Comp2_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_Comp2, "derived_quantities/fitted_Comp2.csv")

#=========================================================================================




#=========================================================================================
#========== Fitted values: Comp3
#=========================================================================================

fitted_Comp3 <- 
  # extract posterior draws for intercept and slope
  as.data.frame(m_diet,
                variable = c("b_Intercept",
                             "b_zi_Intercept",
                             "b_Comp3_z",
                             paste("r_taxon[", taxon_levels, 
                                   ",Comp3_z]", 
                                   sep=""))) %>%
  # assign number to iteration based on row
  mutate(step = row_number()) %>%
  # convert values by taxon to long format
  pivot_longer(4:(length(unique(diet_long_analysis$taxon))+3),
               names_to = "taxon_c",
               values_to = "r_b") %>%
  # add tidy taxon labels
  full_join(
    data.frame(taxon_c = paste("r_taxon[", taxon_levels, 
                               ",Comp3_z]", 
                               sep=""),
               taxon = taxon_levels)) %>%
  # remove verbose taxon labels
  select(-taxon_c) %>%
  # join with taxon-specific zi
  full_join(
    # extract posterior draws for zi by taxon
    as.data.frame(m_diet,
                  variable = paste("r_taxon__zi[", taxon_levels, 
                                   ",Intercept]", 
                                   sep="")) %>%
      # assign number to iteration based on row
      mutate(step = row_number()) %>%
      # convert values by taxon to long format
      pivot_longer(1:length(unique(diet_long_analysis$taxon)),
                   names_to = "taxon_c",
                   values_to = "r_zi") %>%
      # add tidy taxon labels
      full_join(data.frame(taxon_c = paste("r_taxon__zi[", taxon_levels, 
                                           ",Intercept]", 
                                           sep=""),
                           taxon = taxon_levels)) %>%
      # remove verbose taxon labels
      select(-taxon_c)) %>%
  # expand posterior draws for values of predictor
  expand(nesting(taxon,
                 step, 
                 b_Intercept,
                 b_zi_Intercept,
                 b_Comp3_z,
                 r_b,
                 r_zi),
         Comp3_z = seq(min(diet_long_analysis$Comp3_z),
                       max(diet_long_analysis$Comp3_z),
                       length.out = 10)) %>%
  # calculate linear predictor
  # convert inverse logit transform to calculate zi
  # calculated expected value on response scale (zero-inflated Poisson with log-link)
  mutate(lin_pred = b_Intercept + (b_Comp3_z + r_b) * Comp3_z,
         zi = exp(b_zi_Intercept + r_zi) / (1 + exp(b_zi_Intercept + r_zi)),
         pred = exp(lin_pred) * (1 - zi)) %>%
  # group by taxon and levels of predictor for suammry
  group_by(taxon,
           Comp3_z) %>%
  # posterior summary and ungroup
  summarize(fit = median(pred),
            lo = quantile(pred, probs = c(0.16)),
            hi = quantile(pred, probs = c(0.84))) %>%
  ungroup() %>%
  # assign tidy taxon names
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_labels))

# write_csv(fitted_Comp3, "derived_quantities/fitted_Comp3.csv")

#=========================================================================================


