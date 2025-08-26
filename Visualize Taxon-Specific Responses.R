#=========================================================================================
#========== Preliminaries
#=========================================================================================

# Load packages
library(tidyverse)
library(lemon)
library(cowplot)

# import data
diet_long_analysis <- read_csv("data/diet_long_analysis.csv")

# import derived quantities
slope_sd_taxon <- read_csv("derived_quantities/slope_sd_taxon.csv")
intercept_sd <- read_csv("derived_quantities/intercept_sd.csv")
taxon_coefs_long <- read_csv("derived_quantities/taxon_coefs_long.csv")
taxon_interepts_F_relto_M <- read_csv("derived_quantities/taxon_interepts_F_relto_M.csv")
taxon_intercepts_sites <- read_csv("derived_quantities/taxon_intercepts_sites.csv")
taxon_int <- read_csv("derived_quantities/taxon_int.csv")
taxon_zi <- read_csv("derived_quantities/taxon_zi.csv")
fitted_length <- read_csv("derived_quantities/fitted_length.csv")
fitted_gut_length <- read_csv("derived_quantities/fitted_gut_length.csv")
fitted_gap <- read_csv("derived_quantities/fitted_gap.csv")
fitted_raker_length <- read_csv("derived_quantities/fitted_raker_length.csv")
fitted_comp1 <- read_csv("derived_quantities/fitted_comp1.csv")
fitted_comp2 <- read_csv("derived_quantities/fitted_comp2.csv")
fitted_comp3 <- read_csv("derived_quantities/fitted_comp3.csv")

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
taxon_label <- c("Chironomidae (adult)", 
                  "Alona spp.", 
                  "Chrinominae", 
                  "Cyclopideae", 
                  "Mollusca", 
                  "Orthocladinae",     
                  "Ostracoda",          
                  "Cladocerans (other)", 
                  "Stickleback eggs",  
                  "Tanypodinae")

# define plotting data with tidy taxon names
diet_long_plot <- diet_long_analysis %>%
  mutate(taxon_label = factor(taxon,
                               levels = taxon_levels,
                               labels = taxon_label))

# define taxon colors for graphing
taxon_colors <- c("darkgoldenrod3",
                  "darkolivegreen",
                  "darkolivegreen3",
                  "darkseagreen4",
                  "palevioletred3",
                  "deepskyblue1",
                  "deepskyblue4",
                  "dodgerblue3",
                  "lightcoral",
                  "peachpuff3")
names(taxon_colors) <- unique(taxon_coefs_long$taxon_label)

# set theme
theme_set(theme_classic() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       r = 1,
                                       b = 1,
                                       l = 1),
                  legend.margin = margin(t = 0,
                                         r = 0,
                                         b = 0,
                                         l = -4),
                  legend.position = "top",
                  legend.text = element_text(size = 8,
                                             margin = margin(l = 0), 
                                             face = "italic"),
                  legend.key.size = unit(1, "lines"),
                  legend.spacing.x = unit(0.4, "lines"),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.line = element_line(linewidth = 0.25),
                  axis.ticks = element_line(linewidth = 0.25)))

# define breaks for scatter plots
scatter_breaks <- ggh4x::facetted_pos_scales(y = list(
  taxon_label == "Chironomidae (adult)" ~ 
    scale_y_continuous(breaks = c(0, 1, 5, 25), trans = "log1p"),
  taxon_label == "Alona spp." ~ 
    scale_y_continuous(breaks = c(0, 1, 10, 100), trans = "log1p"),  
  taxon_label == "Chrinominae" ~ 
    scale_y_continuous(breaks = c(0, 1, 4, 16), trans = "log1p"),
  taxon_label == "Cyclopideae" ~ 
    scale_y_continuous(breaks = c(0, 1, 10, 100), trans = "log1p"),
  taxon_label == "Mollusca" ~ 
    scale_y_continuous(breaks = c(0, 1, 3, 9), trans = "log1p"),
  taxon_label == "Orthocladinae" ~ 
    scale_y_continuous(breaks = c(0, 1, 6, 36), trans = "log1p"),
  taxon_label == "Ostracoda" ~ 
    scale_y_continuous(breaks = c(0, 1, 3, 9), trans = "log1p"),
  taxon_label == "Cladocerans (other)" ~ 
    scale_y_continuous(breaks = c(0, 1, 6, 36), trans = "log1p"),
  taxon_label == "Stickleback eggs" ~ 
    scale_y_continuous(breaks = c(0, 1, 5, 25), trans = "log1p"),
  taxon_label == "Tanypodinae"~ 
    scale_y_continuous(breaks = c(0, 1, 5, 25), trans = "log1p")
))

#=========================================================================================





#=========================================================================================
#========== SDs characterizing taxon-specific responses
#=========================================================================================

# define factor orders for plotting
slope_sd_taxon_sum <- slope_sd_taxon %>%
  mutate(name = factor(name,
                       levels = c("PC1",
                                  "Gap",
                                  "PC2",
                                  "PC3",
                                  "Raker Length",
                                  "Body Length",
                                  "Gut Length"))) %>%
  group_by(name) %>%
  summarize(Estimate = quantile(value, probs = 0.5),
            Q5 = quantile(value, probs = 0.025),
            Q16 = quantile(value, probs = 0.16),
            Q84 = quantile(value, probs = 0.84),
            Q95 = quantile(value, probs = 0.975))

# define factor orders for plotting
slope_sd_taxon_order <- {slope_sd_taxon_sum %>% arrange(Estimate)}$name
slope_sd_taxon_sum$name <- factor(slope_sd_taxon_sum$name, 
                            levels = slope_sd_taxon_order)

# plot
fig_sds_slope <- ggplot(data = slope_sd_taxon_sum,
                      aes(x = name,
                          y = Estimate))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_point()+
  scale_x_discrete("")+
  scale_y_continuous("Variation among taxon respones (SD)", 
                     expand = c(0, 0),
                     limits = c(0, 1.6),
                     breaks = c(0, 0.5, 1, 1.5))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)))
fig_sds_slope

# export
# ggsave("figures/fig_sds_slope.pdf",
#        plot = fig_sds_slope,
#        width = 3.5,
#        height = 3)

#=========================================================================================




#=========================================================================================
#========== SDs characterizing taxon-specific responses
#=========================================================================================

# summarize posteriors
intercept_sd_sum <- intercept_sd %>%
  mutate(name = factor(name,
                       levels = c("Fish ID",
                                  "Site",
                                  "Taxon",
                                  "Taxon x Site",
                                  "Taxon x Sex"))) %>%
  group_by(name) %>%
  summarize(Estimate = quantile(value, probs = 0.5),
            Q5 = quantile(value, probs = 0.025),
            Q16 = quantile(value, probs = 0.16),
            Q84 = quantile(value, probs = 0.84),
            Q95 = quantile(value, probs = 0.975))

# define factor orders for plotting
intercept_sd_order <- {intercept_sd_sum %>% arrange(Estimate)}$name
intercept_sd_sum$name <- factor(intercept_sd_sum$name, 
                                        levels = intercept_sd_order)

# plot
fig_sds_int <- ggplot(data = intercept_sd_sum,
                      aes(x = name,
                          y = Estimate))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_point()+
  scale_x_discrete("")+
  scale_y_continuous("Variation among intercepts (SD)", 
                     expand = c(0, 0),
                     limits = c(0, 3.1),
                     breaks = c(0, 1, 2, 3))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)))
fig_sds_int

# export
# ggsave("figures/fig_sds_int.pdf",
#        plot = fig_sds_int,
#        width = 3.5,
#        height = 3)




# combine
fig_sds <- plot_grid(NULL,
                     fig_sds_slope,
                     NULL,
                     fig_sds_int,
                     nrow = 4,
                     rel_heights = c(0.075, 1, 0.2, 1),
                     labels = c("","(a)","","(b)"),
                     label_size = 10,
                     label_fontface = "plain",
                     vjust = c(0,-0.4, 0, -0.4),
                     hjust = c(0,-0.05, 0, -0.05))
fig_sds
# ggsave("figures/fig_sds.pdf",
#        plot = fig_sds,
#        width = 3.5,
#        height = 6)

#=========================================================================================





#=========================================================================================
#========== Taxon-specific responses
#=========================================================================================

# panel (a): length responses 
length_responses <- taxon_coefs_long %>% filter(coef == "length_z")
length_responses_order <- {length_responses %>% arrange(Estimate)}$taxon_label
length_responses$taxon_label <- factor(length_responses$taxon_label, 
                                      levels = length_responses_order)
fig_length <- ggplot(data = length_responses,
                   aes(x = taxon_label,
                       y = Estimate,
                       color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("Body length response", 
                     breaks = c(-1, 0, 1),
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_length

# panel (b): gut length responses 
gut_length_responses <- taxon_coefs_long %>% filter(coef == "gut_length_z")
gut_length_responses_order <- {gut_length_responses %>% arrange(Estimate)}$taxon_label
gut_length_responses$taxon_label <- factor(gut_length_responses$taxon_label, 
                                          levels = gut_length_responses_order)
fig_gut_length <- ggplot(data = gut_length_responses,
                           aes(x = taxon_label,
                               y = Estimate,
                               color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("Gut length response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_gut_length

# panel (c): raker length responses 
raker_length_responses <- taxon_coefs_long %>% filter(coef == "raker_length_z")
raker_length_responses_order <- {raker_length_responses %>% arrange(Estimate)}$taxon_label
raker_length_responses$taxon_label <- factor(raker_length_responses$taxon_label, 
                                          levels = raker_length_responses_order)
fig_raker_length <- ggplot(data = raker_length_responses,
                         aes(x = taxon_label,
                             y = Estimate,
                             color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("Raker length response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_raker_length

# panel (d): raker length responses 
comp3_responses <- taxon_coefs_long %>% filter(coef == "Comp3_z")
comp3_responses_order <- {comp3_responses %>% arrange(Estimate)}$taxon_label
comp3_responses$taxon_label <- factor(comp3_responses$taxon_label, 
                                            levels = comp3_responses_order)
fig_comp3 <- ggplot(data = comp3_responses,
                           aes(x = taxon_label,
                               y = Estimate,
                               color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("PC3 response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_comp3

# panel (e): PC2 responses 
comp2_responses <- taxon_coefs_long %>% filter(coef == "Comp2_z")
comp2_responses_order <- {comp2_responses %>% arrange(Estimate)}$taxon_label
comp2_responses$taxon_label <- factor(comp2_responses$taxon_label, 
                                       levels = comp2_responses_order)
fig_comp2 <- ggplot(data = comp2_responses,
                    aes(x = taxon_label,
                        y = Estimate,
                        color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("PC2 response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_comp2

# panel (f): Gap responses 
gap_responses <- taxon_coefs_long %>% filter(coef == "Gap_z")
gap_responses_order <- {gap_responses %>% arrange(Estimate)}$taxon_label
gap_responses$taxon_label <- factor(gap_responses$taxon_label, 
                                       levels = gap_responses_order)
fig_gap <- ggplot(data = gap_responses,
                    aes(x = taxon_label,
                        y = Estimate,
                        color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("Gap response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_gap

# panel (g): PC1 responses 
comp1_responses <- taxon_coefs_long %>% filter(coef == "Comp1_z")
comp1_responses_order <- {comp1_responses %>% arrange(Estimate)}$taxon_label
comp1_responses$taxon_label <- factor(comp1_responses$taxon_label, 
                                       levels = comp1_responses_order)
fig_comp1 <- ggplot(data = comp1_responses,
                    aes(x = taxon_label,
                        y = Estimate,
                        color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("PC1 response", 
                     breaks = c(-1, 0, 1), 
                     limits = c(-1.55, 1.55))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_comp1

# panel (h): Sex differences
sex_responses_order <- {taxon_interepts_F_relto_M %>% arrange(Estimate)}$taxon_label
taxon_interepts_F_relto_M$taxon_label <- factor(taxon_interepts_F_relto_M$taxon_label, 
                                                 levels = sex_responses_order)
fig_sex_taxon <- ggplot(data = taxon_interepts_F_relto_M,
                        aes(x = taxon_label,
                            y = Estimate,
                            color = taxon_label))+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2)+
  geom_point()+
  scale_color_manual("", values = taxon_colors, guide = guide_legend(reverse = TRUE) )+
  scale_x_discrete("")+
  scale_y_continuous("Females vs Males",
                     breaks = c(-3, 0, 3),
                     limits = c(-3.5, 3.5))+
  guides(color = "none")+
  coord_flip()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)),
        axis.text.y = element_text(size = 9))
fig_sex_taxon

# combine
fig_responses <- plot_grid(NULL, NULL, NULL, 
                           fig_gut_length, NULL, fig_length, 
                           NULL, NULL, NULL, 
                           fig_raker_length, NULL, fig_comp3,
                           NULL, NULL, NULL, 
                           fig_comp2, NULL, fig_gap, 
                           NULL, NULL, NULL, 
                           fig_comp1, NULL, fig_sex_taxon,
                           nrow = 8,
                           ncol = 3,
                           rel_heights = c(0.075, 1, 0.15, 1, 0.15, 1, 0.15, 1),
                           rel_widths = c(1, 0.05, 1),
                           labels = c("","","",
                                      "(a)","","(b)",
                                      "","","",
                                      "(c)","","(d)",
                                      "","","",
                                      "(e)","","(f)",
                                      "","","",
                                      "(g)","","(h)"),
                           label_size = 10,
                           label_fontface = "plain",
                           vjust = c(0, 0, 0,
                                     -0.4, 0, -0.4, 
                                     0, 0, 0, 
                                     -0.4,0, -0.4,
                                     0, 0, 0,
                                     -0.4, 0, -0.4, 
                                     0, 0, 0, 
                                     -0.4,0, -0.4),
                           hjust = c(0, 0, 0,
                                     -0.05, 0, -0.05, 
                                     0, 0, 0, 
                                     -0.05,0, -0.05,
                                     0, 0, 0,
                                     -0.05, 0, -0.05, 
                                     0, 0, 0, 
                                     -0.05,0, -0.05))
fig_responses

# export
ggsave("figures/fig_responses.pdf",
       plot = fig_responses,
       width = 6.5,
       height = 9)

#=========================================================================================






#=========================================================================================
#========== Spatial variation in taxonomic representation in diet
#=========================================================================================

# extract taxon variation among sites
taxon_sites <- tibble(rownames_to_column(
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
taxon_sites_clean <- taxon_sites %>% 
  mutate(taxon = rowname %>% map_chr(~taxon_f(.x)),
         Site = rowname %>% map_chr(~site_f(.x))) %>%
  select(taxon, Site, Estimate, Est.Error, Q5, Q95, Q16, Q84) %>%
  mutate(taxon_label = factor(taxon,
                               levels = c("adult_fly", 
                                          "alona", 
                                          "chrinominae", 
                                          "copepods", 
                                          "mollusks", 
                                          "orthocladinae",     
                                          "ostracod",          
                                          "other_cladocerans", 
                                          "stickleback_eggs",  
                                          "tanypodinae"),
                               labels = c("Chironomidae (adult)", 
                                          "Alona spp.", 
                                          "Chrinominae", 
                                          "Cyclopideae", 
                                          "Mollusca", 
                                          "Orthocladinae",     
                                          "Ostracoda",          
                                          "Cladocerans (other)", 
                                          "Stickleback Eggs",  
                                          "Tanypodinae")))


taxon_intercepts_sites <- taxon_intercepts_sites %>%
  mutate(taxon_label = factor(taxon_label,
                               levels = taxon_labels))

taxon_int <- taxon_int %>%
  mutate(taxon_label = factor(taxon_label,
                               levels = taxon_labels))

# graph
ggplot(data = taxon_intercepts_sites,
       aes(x = Site,
           y = Estimate,
           color = taxon_label))+
  facet_wrap(~taxon_label, scales = "free_y", ncol = 2)+
  geom_errorbar(aes(ymin = Q16,
                    ymax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(ymin = Q5,
                    ymax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_hline(data = taxon_int,
             aes(yintercept = Estimate), 
             linetype = 2, 
             linewidth = 0.2)+
  geom_point()+
  guides(color = "none")+
  scale_color_manual("", values = taxon_colors)+
  coord_flip()+
  scale_x_discrete("")+
  scale_y_continuous("Intercept estimate", 
                     breaks = c(-5, -2, 1))+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

# ggsave("figures/fig_stations.pdf", 
#        plot = last_plot(), 
#        width = 5,
#        height = 7)

#=========================================================================================
