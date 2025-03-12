
library(RInSp)
library(vegan)
library(MCMCglmm)
library(tidyverse)
library(brms)

setwd("/Users/kashastrickland/Desktop/diet_specialisation")

############Joes data cleaning steps

# Import data
diet <- read.csv("Final.csv")

# Gather data in long format
# Select key variables
# Fill blanks with 0's
diet_long <- diet %>% 
  select(fish_id,
         length,
         gut_length,
         Gap,
         NR.2,
         NR.3,
         Station,
         basin,
         tanypodinae,
         orthocladinae,
         chironominii,
         tanytarsini,
         unknown.larvae,
         adult_fly,
         Alona.rectangula,
         Alona.affinis.quad,
         Unknown.alona,
         alonella.nana,
         macrothrix.hirsuticornix,
         acroperus.harpae,
         daphnia.longispina,
         resting.eggs,
         Copepods_adult,
         Copepods_nauplii,
         ostracod,
         snail,
         bivalve,
         stickleback.eggs
  ) %>%
  pivot_longer(cols = c(tanypodinae,
                        orthocladinae,
                        chironominii,
                        tanytarsini,
                        unknown.larvae,
                        adult_fly,
                        Alona.rectangula,
                        Alona.affinis.quad,
                        Unknown.alona,
                        alonella.nana,
                        macrothrix.hirsuticornix,
                        acroperus.harpae,
                        daphnia.longispina,
                        resting.eggs,
                        Copepods_adult,
                        Copepods_nauplii,
                        ostracod,
                        snail,
                        bivalve,
                        stickleback.eggs),
               names_to = "taxon",
               values_to = "number") %>%
  mutate(number = ifelse(is.na(number), 0, number)) %>%
  na.omit()

# For the taxon specific analysis, we need to select taxa that have 
#   enough observations / aggregate to get to that level
# To facilitate this, it is worth looking at the proportion of fish 
#   in which each taxon is found
diet_long %>%
  mutate(present = ifelse(number > 0, 1, 0)) %>%
  group_by(taxon) %>%
  summarize(mean = mean(present))

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
         alona = Alona.affinis.quad  + Alona.rectangula + Unknown.alona,
         copepods = Copepods_adult + Copepods_nauplii,
         other_cladocerans = acroperus.harpae + alonella.nana + daphnia.longispina + 
           macrothrix.hirsuticornix,
         mollusks = snail + bivalve,
         Site = ifelse(Station == "Lake" & basin == "North", "Lake/North",
                       ifelse(Station == "Lake" & basin == "South", "Lake/South",
                              ifelse(Station == "Shore", "Shore/North", Station)))) %>%
  rename(stickleback_eggs = stickleback.eggs)

######
##variation in diet across sites

##tidy up data
diet_clean<-subset(diet_wide, !rowSums(diet_wide[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                                                    "alona","other_cladocerans","stickleback_eggs","copepods",
                                                    "mollusks","ostracod")])==0)

dietR<-diet_clean[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                  "alona","other_cladocerans","stickleback_eggs","copepods",
                  "mollusks","ostracod")]

names(dietR)<-c("Tanypodinae","Orthocladinae","Chrinominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Cyclopideae","Mollusca","Ostracoda")

colnames(diet_clean)[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                       "alona","other_cladocerans","stickleback_eggs","copepods",
                       "mollusks","ostracod")]<-c("Tanypodinae","Orthocladinae","Chrinominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Cyclopideae","Mollusca","Ostracoda")

diet_clean$n_items<-rowSums(dietR)

###Run RDA ordination
st_rda<-rda(dietR~
      Site+length,data=diet_clean,scale=TRUE)

fit <- envfit(st_rda, dietR, perm = 999)

col_palette <- c("#f0e442","#0072b2","#d55e00", "#cc79a7")
col_palette

pca_scores<-scores(st_rda)
pdf("Figure4.pdf",width=11,height=9)
plot(pca_scores$sites[,1],
     pca_scores$sites[,2],
     pch=21,
     bg=as.numeric(diet_clean$Site),
     xlim=c(-3,5.5),
     xlab="RDA1",
     ylab="RDA2")
plot(fit, col = "red")
ordiellipse(st_rda,diet_clean$Site,col = col_palette,conf=0.8)
legend('topright', legend=unique(diet_clean$Site), col=unique(col_palette), pch = 16)
dev.off()

##################
###individual specialisation across the lake
lkWide_spec<-import.RInSp(dietR)
lake_niche<-WTcMC(lkWide_spec,weight="N_items",replicates=1500)
psi_lkwide<-PSicalc(lkWide_spec,pop.diet = "average",replicates=1500)
hist(psi_lkwide$PSi,xlab="PSi")
diet_clean$PSi<-psi_lkwide$PSi

##specialisation according to location
###relationship between PSi and morphology
# gut length
m_gut <- lm(gut_length ~ length, data = diet_clean)
summary(m_gut)
diet_clean$gut_length_res <- resid(m_gut)
# Gap 
m_gap <- lm(Gap ~ length, data = diet_clean)
summary(m_gap)
diet_clean$Gap_res <- resid(m_gap)
# Gill raker 2
m_raker <- lm(NR.2 ~ length, data = diet_clean)
summary(m_raker)
diet_clean$NR.2_res <- resid(m_raker)

str(diet_clean)
diet_clean <- diet_clean %>%
  mutate(length_z = c(scale(length)),
         gut_length_z = c(scale(gut_length)),
         Gap_z = c(scale(Gap)),
         NR.2_z = c(scale(NR.2)))

m1<-lm(PSi~length_z+gut_length_z+Gap_z+NR.2_z+Site,
              data=diet_clean)
summary(m1)
plot(m1)

pdf("Figure5a.pdf",width=9)
ggplot(diet_clean,aes(y=PSi,x=Site,fill=Site))+
  geom_boxplot()+
  ggsignif::geom_signif(comparisons=list(c("HS","Lake/South"),c("HS","Shore/North")),step_increase = 0.1,  annotations ="*") +
  xlab("\nSite")+
  ylab("PSi\n") +
  scale_fill_manual(values = c("#d55e00", "#cc79a7", "#0072b2", "#f0e442"),guide="none")+
  theme_classic()+
  theme(text=element_text(size=16))
dev.off()

##include head shape PCs
list.files()
pc_hs<-read.csv("pca_scores.csv")
colnames(pc_hs)[1]<-"fish_id"
diet_clean_hs<-left_join(diet_clean,pc_hs[1:4])

str(diet_clean_hs)
pairs(diet_clean_hs[c("length_z","gut_length_z","Gap_z","NR.2_z","Comp1","Comp2","Comp3")])
##add sex in for model
diet_clean_hs<-left_join(diet_clean_hs,diet[,c("fish_id","sex")])
##run model to explore differences in PSi from morphology or site
library(brms)
m2<-brm(PSi~length_z+gut_length_z+Gap_z+NR.2_z+Comp1+Comp2+Comp3+sex+(1|Site),
        chains = 4,
        iter = 4000,
        warmup = 2000,
        control =list(adapt_delta=0.99), 
        data=diet_clean_hs)

summary(m2)
plot(m2)

###is specialisation linked to certain diets and morphology?
diet_psi<-diet_clean[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                       "alona","other_cladocerans","stickleback_eggs","copepods",
                       "mollusks","ostracod","PSi")]

names(diet_psi)<-c("Tanypodinae","Orthocladinae","Chrinominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Cyclopideae","Mollusca","Ostracoda","PSi")
st_rda_psi<-rda(diet_psi~
              Site+length,data=diet_clean,scale=TRUE)
plot(st_rda_psi)
summary(st_rda_psi)
anova(st_rda_psi,by="term", permutations=199)
fit_psi <- envfit(st_rda_psi, diet_psi, perm = 999, choices = c(1,2))
plot(st_rda_psi, type = "p",scaling=2)
plot(fit_psi)
pca_scores_psi<-scores(st_rda_psi)
pdf("Figure5b.pdf",width=12,height=9)
plot(pca_scores_psi$sites[,1],
     pca_scores_psi$sites[,2],
     pch=21,
     bg=as.numeric(diet_clean$Site),
     xlim=c(-3,5.5),
     xlab="RDA1",
     ylab="RDA2")
#ylim=c(-2,2))
plot(fit_psi, col = "red")
ordiellipse(st_rda_psi,diet_clean$Site,col = col_palette,conf=0.8)
legend('topright', legend=unique(diet_clean$Site), col=unique(col_palette), pch = 16)
dev.off()

###HOT SHORE
dietHS<-diet_clean[diet_clean$Site=="HS",]

###individual specialisation in HS
HS_spec<-import.RInSp(dietHS[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                                   "alona","other_cladocerans","stickleback_eggs","copepods",
                                   "mollusks","ostracod")])
HS_niche<-WTcMC(HS_spec,weight="N_items",replicates=1500)
psi_HS<-PSicalc(HS_spec,pop.diet = "average",replicates=1500)
hist(psi_HS$PSi)
dietHS$PSiHS<-psi_HS$PSi

###Shore North
dietNS<-diet_clean[diet_clean$Site=="Shore/North",]
###individual specialisation in NS
NS_spec<-import.RInSp(dietNS[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                               "alona","other_cladocerans","stickleback_eggs","copepods",
                               "mollusks","ostracod")])
NS_niche<-WTcMC(NS_spec,weight="N_items",replicates=1500)
psi_NS<-PSicalc(NS_spec,pop.diet = "average",replicates=1500)
hist(psi_NS$PSi)
dietNS$PSiNS<-psi_NS$PSi

###Lake North
dietNL<-diet_clean[diet_clean$Site=="Lake/North",]

###individual specialisation in NL
NL_spec<-import.RInSp(dietNL[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                               "alona","other_cladocerans","stickleback_eggs","copepods",
                               "mollusks","ostracod")])

NL_niche<-WTcMC(NL_spec,weight="N_items",replicates=1500)

psi_NL<-PSicalc(NL_spec,pop.diet = "average",replicates=1500)
hist(psi_NL$PSi)
dietNL$PSiNL<-psi_NL$PSi

###Lake South
dietSL<-diet_clean[diet_clean$Site=="Lake/South",]

###individual specialisation in HS
SL_spec<-import.RInSp(dietSL[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                               "alona","other_cladocerans","stickleback_eggs","copepods",
                               "mollusks","ostracod")])

SL_niche<-WTcMC(SL_spec,weight="N_items",replicates=1500)
psi_SL<-PSicalc(SL_spec,pop.diet = "average",replicates=1500)
hist(psi_SL$PSi)
dietSL$PSiSL<-psi_SL$PSi

##morphological variation across sites
morph<-left_join(diet,pc_hs[1:4])
str(morph)  
morph$Site<-ifelse(morph$Station == "Lake" & morph$basin == "North", "Lake/North",
                          ifelse(morph$Station == "Lake" & morph$basin == "South", "Lake/South",
                                 ifelse(morph$Station == "Shore", "Shore/North", morph$Station)))


mvl<-bf(scale(length)~sex+(1|p|Site))
mvgl<-bf(scale(gut_length)~sex+scale(length)+(1|p|Site))
mvgr2<-bf(scale(NR.2)~sex+scale(length)+(1|p|Site))
mvgap<-bf(scale(Gap)~sex+scale(length)+(1|p|Site))
mvPC1<-bf(scale(Comp1)~sex+scale(length)+(1|p|Site))
mvPC2<-bf(scale(Comp2)~sex+scale(length)+(1|p|Site))
mvPC3<-bf(scale(Comp3)~sex+scale(length)+(1|p|Site))

fit1 <- brm(mvl+mvgl+mvgr2+mvgap+mvPC1+mvPC2+mvPC3+set_rescor(TRUE), 
            data = morph, control =list(adapt_delta=0.99,max_treedepth=15),chains = 4, cores = 4)
summary(fit1)
plot(fit1)

##extract site BLUPs with confidence intervals and plot
rfLength<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 1])[c(1,3:6)]
rfgutLength<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 2])[c(1,3:6)]
rfGR2<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 3])[c(1,3:6)]
rfGap<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 4])[c(1,3:6)]
rfComp1<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 5])[c(1,3:6)]
rfComp2<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 6])[c(1,3:6)]
rfComp3<-as.data.frame(coef(fit1,probs = c(0.05, 0.16, 0.84, 0.95))$Site[ , , 7])[c(1,3:6)]

rfLength$site<-rownames(rfLength)
rfgutLength$site<-rownames(rfgutLength)
rfGR2$site<-rownames(rfGR2)
rfGap$site<-rownames(rfGap)
rfComp1$site<-rownames(rfComp1)
rfComp2$site<-rownames(rfComp2)
rfComp3$site<-rownames(rfComp3)

rfLength$trait<-"Body length"
rfgutLength$trait<-"Gut length"
rfGR2$trait<-"Gill raker length"
rfGap$trait<-"Gill raker gap width" 
rfComp1$trait<-"PC1"
rfComp2$trait<-"PC2"
rfComp3$trait<-"PC3"

rf_all<-rbind(rfLength,rfgutLength,rfGR2,rfGap,rfComp1,rfComp2,rfComp3)

##get trait specific intercepts
int_traits<-data.frame(trait=unique(rf_all$trait),estimate=fixef(fit1)[1:7])

# set theme
theme_set(theme_classic() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_text(size=10),
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

col_paletteM <- c("#d55e00", "#cc79a7","#0072b2","#f0e442")

pdf("Figure1.pdf")
ggplot(rf_all,aes(x=Estimate, y = site,color=site))  +
  geom_point()+
  facet_wrap(~trait, scales = "free", ncol = 2)+
  geom_errorbar(aes(xmin = Q16,
                    xmax = Q84),
                width = 0,
                linewidth = 0.75)+
  geom_errorbar(aes(xmin = Q5,
                    xmax = Q95),
                width = 0,
                linewidth = 0.25)+
  geom_vline(data = int_traits,
             aes(xintercept = estimate), 
             linetype = 2, 
             linewidth = 0.2)+
  guides(color = "none") +
  scale_x_continuous("Intercept estimate (z-scale) \n",limits=c(-1,1))+
  scale_y_discrete("")+
  scale_color_manual(values=col_paletteM)
dev.off()

