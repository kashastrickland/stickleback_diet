#############################
## Author: Kasha Strickland, 2025

library(RInSp)
library(vegan)
library(MCMCglmm)
library(tidyverse)
library(brms)

setwd("")

######
diet_wide<-read_csv("data/diet_wide.csv")
pc_hs<-read.csv("data/pca_scores.csv")

#############################################
##morphological variation across sites
morph<-left_join(diet_wide,pc_hs[1:4])
str(morph)  

mvl<-bf(length_z~sex+(1|p|Site))
mvgl<-bf(gut_length_z~sex+length_z+(1|p|Site))
mvgr2<-bf(raker_length_z~sex+length_z+(1|p|Site))
mvgap<-bf(Gap_z~sex+length_z+(1|p|Site))
mvPC1<-bf(Comp1_z~sex+length_z+(1|p|Site))
mvPC2<-bf(Comp2_z~sex+length_z+(1|p|Site))
mvPC3<-bf(Comp3_z~sex+length_z+(1|p|Site))

fit1 <- brm(mvl+mvgl+mvgr2+mvgap+mvPC1+mvPC2+mvPC3+set_rescor(TRUE), 
            data = morph, control =list(adapt_delta=0.99,max_treedepth=15),chains = 4, cores = 4)
summary(fit1)
plot(fit1)

##make results table
morph_fe<-brms::fixef(fit1,probs = c(0.05,0.95))[,c(1,3:4)]
morph_fe<-data.frame(parameter=rownames(morph_fe),round(morph_fe,3))
morph_fe$trait<-vapply(strsplit(morph_fe$parameter,"_"), `[`, 1, FUN.VALUE=character(1))
morph_fe$parameter<-vapply(strsplit(morph_fe$parameter,"_"), `[`, 2, FUN.VALUE=character(1))
str(morph_fe)
morph_feSPL<-split(morph_fe,morph_fe$parameter)
morph_feSPL[[3]]
table1<-data.frame(Trait=c("Body length","Gut length","Gill raker length","Gill raker gap width", "PC1","PC2","PC3"),
                    Intercept=paste0(morph_feSPL[[1]]$Estimate," (",morph_feSPL[[1]]$Q5,"; ",morph_feSPL[[1]]$Q95,")"),
                    SexM=paste0(morph_feSPL[[3]]$Estimate," (",morph_feSPL[[3]]$Q5,"; ",morph_feSPL[[3]]$Q95,")"),
                    Bodylength=paste0(c("",morph_feSPL[[2]]$Estimate)," (",c("",morph_feSPL[[2]]$Q5),"; ",c("",morph_feSPL[[2]]$Q95),")"))

morph_re<-brms::VarCorr(fit1,probs = c(0.05,0.95))
table1$Site<-paste0(round(morph_re$Site$sd[,1],3)," (",round(morph_re$Site$sd[,2],3),"; ",round(morph_re$Site$sd[,3],3),")")
table1$REsidual<-paste0(round(morph_re$residual__$sd[,1],3)," (",round(morph_re$residual__$sd[,2],3),"; ",round(morph_re$residual__$sd[,3],3),")")
write.table(table1,"tables_figures/table1.txt",sep="\t",row.names = F,quote=F)


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

col_paletteM <- c("#d55e00", "#cc79a7","#0072b2","darkgoldenrod1")

pdf("tables_figures/Figure3.pdf")
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

#####################
##tidy up data by removing individuals with empty stomachs
diet_clean<-subset(diet_wide, !rowSums(diet_wide[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                                                   "alona","other_cladocerans","stickleback_eggs","copepods",
                                                   "mollusks","ostracod")])==0)

dietR<-diet_clean[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                  "alona","other_cladocerans","stickleback_eggs","copepods",
                  "mollusks","ostracod")]

names(dietR)<-c("Tanypodinae","Orthocladinae","Chironominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Copepoda","Mollusca","Ostracoda")

colnames(diet_clean)[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                       "alona","other_cladocerans","stickleback_eggs","copepods",
                       "mollusks","ostracod")]<-c("Tanypodinae","Orthocladinae","Chironominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Copepoda","Mollusca","Ostracoda")

diet_clean$n_items<-rowSums(dietR)

##diet_clean = empty stomachs removed

###Run RDA ordination
st_rda<-rda(dietR~
      Site+length,data=diet_clean,scale=TRUE)

fit <- envfit(st_rda, dietR, perm = 999)

col_palette <- c("darkgoldenrod1","#0072b2","#d55e00", "#cc79a7")
col_palette

pca_scores<-scores(st_rda)
pdf("tables_figures/FigureS1.pdf",width=14,height=13)
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

pdf("tables_figures/Figure6a.pdf",width=9)
ggplot(diet_clean,aes(y=PSi,x=Site,fill=Site))+
  geom_boxplot()+
  ggsignif::geom_signif(comparisons=list(c("HS","Lake/South"),c("HS","Shore/North")),step_increase = 0.1,  annotations ="*") +
  xlab("\nSite")+
  ylab("PSi\n") +
  scale_x_discrete(labels=c("Warm springs", "North basin lake","South basin lake","North basin shore"))+
  scale_fill_manual(values = c("#d55e00", "#cc79a7", "#0072b2", "darkgoldenrod1"),guide="none")+
  theme_classic()+
  theme(text=element_text(size=16))
dev.off()

###run model to explore differences in PSi from morphology or site
m1<-brm(PSi~length_z+gut_length_z+raker_length_z+Gap_z+Comp1+Comp2+Comp3+sex+(1|Site),
        chains = 4,
        iter = 6000,
        warmup = 3000,
        control =list(adapt_delta=0.9999,max_treedepth=15), 
        data=diet_clean)

summary(m1)
plot(m1)

##make a results table
psi_fe<-brms::fixef(m1,probs = c(0.05,0.95))[,c(1,3:4)]
psi_fe<-data.frame(parameter=rownames(psi_fe),round(psi_fe,3))
psi_re<-brms::VarCorr(m1,probs = c(0.05,0.95))
psi_re<-data.frame(parameter=c("Site","Residual"),Estimate=c(psi_re[[1]]$sd[,1],psi_re[[2]]$sd[,1]),
                   Q5=c(psi_re[[1]]$sd[,3],psi_re[[2]]$sd[,3]),Q95=c(psi_re[[1]]$sd[,4],psi_re[[2]]$sd[,4]))

psi_re[2:4]<-round(psi_re[2:4],3)

table4<-rbind(psi_fe,psi_re)
table4<-data.frame(Parameter=table4$parameter,
                   Estimate = paste0 (table4$Estimate," (",table4$Q5,"; ",table4$Q95,")"))
write.table(table4,"tables_figures/table4.txt",sep="\t",row.names = F,quote=F)

###is specialisation linked to certain diets and morphology?
diet_psi<-diet_clean[c("tanypodinae","orthocladinae","chironominii","adult_fly",
                       "alona","other_cladocerans","stickleback_eggs","copepods",
                       "mollusks","ostracod","PSi")]

names(diet_psi)<-c("Tanypodinae","Orthocladinae","Chironominae","Chironomidae (adult)","Alona spp.","Cladocerans (other)","Stickleback eggs","Copepoda","Mollusca","Ostracoda","PSi")
st_rda_psi<-rda(diet_psi~
              Site+length,data=diet_clean,scale=TRUE)
plot(st_rda_psi)
summary(st_rda_psi)
anova(st_rda_psi,by="term", permutations=199)
fit_psi <- envfit(st_rda_psi, diet_psi, perm = 999, choices = c(1,2))
plot(st_rda_psi, type = "p",scaling=2)
plot(fit_psi)
pca_scores_psi<-scores(st_rda_psi)
pdf("tables_figures/Figure6b.pdf",width=12,height=9)
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
