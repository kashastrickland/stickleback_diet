#=========================================================================================
#========== Preliminaries
#=========================================================================================

# Load packages
library(tidyverse)

#=========================================================================================





#=========================================================================================
#========== Data Preparation
#=========================================================================================

diet_import <- read_csv("data/dissection_data.csv")
pca_scores <- read_csv("data/pca_scores.csv")

# add first three pc axes to diet data
diet <- full_join(diet_import, pca_scores %>% select(fish_id, Comp1, Comp2, Comp3))


# Gather data in long format
# Average gill raker length for 2nd and 3rd rakers
# Select key variables
# Fill blanks with 0's
diet_long <- diet %>% 
  mutate(raker_length = (NR.2 + NR.3) / 2) %>%
  select(fish_id,
         sex,
         length,
         gut_length,
         Gap,
         raker_length,
         Comp1,
         Comp2,
         Comp3,
         Station,
         basin,
         tanypodinae,
         orthocladinae,
         chironominii,
         tanytarsini,
         `unknown larvae`,
         adult_fly,
         `Alona rectangula`,
         `Alona affinis/quad`,
         `Unknown alona`,
         `alonella nana`,
         `macrothrix hirsuticornix`,
         `acroperus harpae`,
         `daphnia longispina`,
         `resting eggs`,
         `Copepods_adult`,
         `Copepods_nauplii`,
         `ostracod`,
         `snail`,
         `bivalve`,
         `stickleback eggs`
  ) %>%
  pivot_longer(cols = c(tanypodinae,
                        orthocladinae,
                        chironominii,
                        tanytarsini,
                        `unknown larvae`,
                        adult_fly,
                        `Alona rectangula`,
                        `Alona affinis/quad`,
                        `Unknown alona`,
                        `alonella nana`,
                        `macrothrix hirsuticornix`,
                        `acroperus harpae`,
                        `daphnia longispina`,
                        `resting eggs`,
                        `Copepods_adult`,
                        `Copepods_nauplii`,
                        `ostracod`,
                        `snail`,
                        `bivalve`,
                        `stickleback eggs`),
               names_to = "taxon",
               values_to = "number") %>%
  mutate(number = ifelse(is.na(number), 0, number)) %>%
  na.omit()

# For the taxon specific analysis, we need to select taxa that have 
#   enough observations / aggregate to get to that level
# To facilitate this, it is worth looking at the proportion of fish 
#   in which each taxon is found
observed_prevalence <- diet_long %>%
  mutate(present = ifelse(number > 0, 1, 0)) %>%
  group_by(taxon) %>%
  summarize(prevalence = mean(present),
            abundance = mean(number),
            sd = sd(number)) %>%
  arrange(-prevalence)
# write_csv(observed_prevalence, "diet_prevalence.csv")


# Three of the major midge groups (Tanypodinae, Orthocladinae, and Chironomini)
#   plus adult midges meet that threshold
# However, the Tanytarsini do not. 
# So, I decided to combine the Chironomini and Tanytarsini into a Chironominae group
# Unknown midge larvae were only found in 2% of fish, and it is hard to know how to combine them with the others
# So, I think they should be left out. 

# Alona rectangula is present in almost 7%.
# But, I think it makes sense to group the Alona's together since 
#   since unknown alonas were present in 20% of fish

# None of the other cladocerans were found in 10% of samples, so I decided to combine them

# Snails were found in nearly 7% of fish. Only 3% of fish had bivalves.
# I decided to combine these as "mollusks" 

# Copepod adults were found in 20% of samples, while the nauplii were only found in 2%
# So, I decided to combine these.

# Finally, I decided to create a "site variable
# This is essentially the same as station, but separating lake stations for north and south 

# prepare wide data frame
diet_wide <- diet_long %>%
  pivot_wider(names_from = taxon,
              values_from = number) %>%
  na.omit() %>%
  mutate(chrinominae = chironominii + tanytarsini,
         alona = `Alona affinis/quad`  + `Alona rectangula` + `Unknown alona`,
         copepods = Copepods_adult + Copepods_nauplii,
         other_cladocerans = `acroperus harpae` +`alonella nana` + `daphnia longispina` + 
           `macrothrix hirsuticornix`,
         mollusks = snail + bivalve,
         Site = ifelse(Station == "Lake" & basin == "North", "Lake/North",
                       ifelse(Station == "Lake" & basin == "South", "Lake/South",
                              ifelse(Station == "Shore", "Shore/North", Station)))) %>%
  rename(stickleback_eggs = `stickleback eggs`)

# When examining relationships with morphological traits, colinearity is a potential concern 

# Check prevalance of chironomids
# diet_wide %>%
#   mutate(Chiro = chrinominae + orthocladinae + adult_fly + `unknown larvae` + tanypodinae,
#          present = ifelse(Chiro > 0, 1, 0)) %>%
#   select(fish_id, Chiro, present) %>%
#   summarize(Chiro_prev = mean(present))

# Define relative length of gut, Gap, and gill rakers as residuals from regression with length
# Focus on gill raker 2 (since it is strongly correlated with gill raker 3)

# gut length
m_gut <- lm(gut_length ~ length, data = diet_wide)
summary(m_gut)
diet_wide$gut_length_res <- resid(m_gut)

# Gap
m_gap <- lm(Gap ~ length, data = diet_wide)
summary(m_gap)
diet_wide$Gap_res <- resid(m_gap)

# Gill raker 
m_raker <- lm(raker_length ~ length, data = diet_wide)
summary(m_raker)
diet_wide$raker_length_res <- resid(m_raker)

# Z-score predictors
diet_wide <- diet_wide %>%
  mutate(length_z = c(scale(length)),
         gut_length_z = c(scale(gut_length_res)),
         Gap_z = c(scale(Gap_res)),
         raker_length_z = c(scale(raker_length_res)),
         Comp1_z = c(scale(Comp1)),
         Comp2_z = c(scale(Comp2)),
         Comp3_z = c(scale(Comp3)))

# examine trait correaltions (to avoid collinearity issues)
round(cor(diet_wide %>% select(length_z, gut_length_z, Gap_z, raker_length_z,
                               Comp1_z, Comp2_z, Comp3_z)),2)

# convert back to long form for mixed models
diet_long_analysis <- diet_wide %>%
  select(fish_id, sex, length, length_z, gut_length, gut_length_z, Gap, Gap_z, 
         raker_length, raker_length_z, Comp1_z, Comp2_z, Comp3_z, Station, basin, Site, 
         tanypodinae,
         orthocladinae,
         chrinominae,
         adult_fly,
         alona,
         other_cladocerans,
         stickleback_eggs,
         copepods,
         mollusks,
         ostracod) %>%
  pivot_longer(cols = c(tanypodinae,
                        orthocladinae,
                        chrinominae,
                        adult_fly,
                        alona,
                        other_cladocerans,
                        stickleback_eggs,
                        copepods,
                        mollusks,
                        ostracod),
               names_to = "taxon",
               values_to = "count_in_stomach") %>%
  select(fish_id,
         sex,
         Station,
         basin,
         Site,
         length,
         length_z,
         gut_length,
         gut_length_z,
         Gap,
         Gap_z,
         raker_length,
         raker_length_z,
         Comp1_z,
         Comp2_z,
         Comp3_z,
         taxon,
         count_in_stomach)

# export data for analysis
write_csv(diet_long_analysis, "data/diet_long_analysis.csv")
write_csv(diet_wide, "data/diet_wide.csv")

#=========================================================================================

##create a figure and table to describe overall prevalence and abundance of diet items across the population
str(diet_long_analysis)

observed_prevalence_concat <- diet_long_analysis %>%
  mutate(present = ifelse(count_in_stomach > 0, 1, 0)) %>%
  group_by(taxon) %>%
  summarize(prevalence = mean(present),
            abundance = mean(count_in_stomach)) %>%
  arrange(+prevalence)

##stacked bar plot

bp_ids<-diet_long_analysis[c("fish_id","taxon","count_in_stomach","Site")]

bp_ids<-bp_ids  %>%
  group_by(fish_id)  %>%
  mutate(percent = count_in_stomach/sum(count_in_stomach)) %>%
  arrange(desc(percent)) %>%
  mutate(sample = forcats::fct_inorder(fish_id))

##remove empty stomachs  
bp_ids <- na.omit(bp_ids)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "copepods", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "mollusks", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "alona", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "chrinominae", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "orthocladinae", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "adult_fly", 0, bp_ids$percent), FUN = sum)
bp_ids$fish_id <-reorder(bp_ids$fish_id,
                         ifelse(!bp_ids$taxon %in% "stickleback_eggs", 0, bp_ids$percent), FUN = sum)

taxon_colors <- c("deepskyblue4",
                  "dodgerblue3",
                  "palevioletred3",
                  "lightcoral",
                  "peachpuff3",
                  "darkseagreen4",
                  "darkolivegreen",
                  "darkolivegreen3",
                  "deepskyblue1",
                  "darkgoldenrod3")

pdf("tables_figures/Fig1B_dietcomp.pdf",width=12,height=9)
ggplot(bp_ids, aes(fill=factor(taxon,levels=observed_prevalence_concat$taxon), y=percent, x= fish_id)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, ncol = 2, scales = "free", labeller = as_labeller(c("HS" = "Warm springs", "Lake/North" = "North basin lake" , "Lake/South" = "South basin lake", "Shore/North"="North basin shore"))) +
  scale_fill_manual(name="",values=taxon_colors,
                    labels=c("Ostracoda","Cladocerans (other)","Mollusca","Stickleback eggs","Tanypodinae",
                             "Copepoda","Alona spp.","Chironominae","Orthocladinae","Chironomidae (adult)"))+
  xlab("")+
  ylab("Relative abundance\n")+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),
        axis.title.y =element_text(color = "grey20", size = 14, face = "plain"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "grey20", size = 14, hjust = .5, vjust = .5, face = "plain"),
        legend.text = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        strip.background = element_blank(),
        strip.text = element_text(color = "grey20", size = 14, face = "plain"))
dev.off()
