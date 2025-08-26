### Set working director

### Load Libraries
#install.packages("geomorph")
library(geomorph)
library(abind) # Combines multidimensional data frames like my tps files
library(ggpubr)
library(ggplot2)
library(plotrix)

### Load all tps files
### These files are raw data from tpsDig, are bent
tps_1 <- readland.tps("data/Loka_GS.TPS", specID = "ID") # Site GS
tps_2 <- readland.tps("data/Loka_NS.TPS", specID = "ID") # Site NS
tps_3 <- readland.tps("data/Loka_4423.TPS", specID = "ID") # Site 44 and 23
tps_4 <- readland.tps("data/Loka_HS2.TPS", specID = "ID") # Site HS
tps_5 <- readland.tps("data/Loka_135.TPS", specID = "ID") # Site 135
tps_6 <- readland.tps("data/Loka_124.TPS", specID = "ID") # Site 124

all <- abind(tps_1,tps_2,tps_3,tps_4,tps_5,tps_6) # Combine all six
dim(all)

## Remove landmarks 5-18 as the fish were bent, only looking at the head
all.new <- all[-5,,] #Removing landmark 5
all.new <- all.new[-5,,] #Removing landmark 6
all.new <- all.new[-5,,] #Removing landmark 7
all.new <- all.new[-5,,] #Removing landmark 8
all.new <- all.new[-5,,] #Removing landmark 9
all.new <- all.new[-5,,] #Removing landmark 10
all.new <- all.new[-5,,] #Removing landmark 11
all.new <- all.new[-5,,] #Removing landmark 12
all.new <- all.new[-5,,] #Removing landmark 13
all.new <- all.new[-5,,] #Removing landmark 14
all.new <- all.new[-5,,] #Removing landmark 15
all.new <- all.new[-5,,] #Removing landmark 16
all.new <- all.new[-5,,] #Removing landmark 17
all.new <- all.new[-5,,] #Removing landmark 18
dim(all.new)

### Load and prep meta data
# The file with all fish with landmarks
info <- read.csv("data/Landmarked.221.csv", sep=";", stringsAsFactors = FALSE, header=TRUE, dec=",")

### Get list of fish id's with open mouth
open_mouth <- subset(info, mouth == "open")$fish_id
### Remove space at the end
### Can't remove and estimate missing landmarks if there is space
dimnames(all.new)[[3]] <- gsub(" ", "", dimnames(all.new)[[3]])

### Start with landmark stuff
### Remove landmarks 1 and 28 for fish with open mouth
all.new[c(1,14),,which(dimnames(all.new)[[3]] %in% open_mouth)] <- NA

### Estimate missing landmarks
all.new <- estimate.missing(all.new, method = "TPS") 

### Define sliding landmarks and save them
### I had done the sliding landmarks all ready so I had the csv ready
#sliders_head <- define.sliders(all.new, nsliders=23)
sliders_head <- read.csv("data/curveslide.csv")

gpa <- gpagen(all.new, curve = sliders_head)
### Find the outliers
outliers <- plotOutliers(gpa$coords) 
### The first five are the outliers
outliers

pc <- gm.prcomp(gpa$coords)
plot(pc)
summary(pc) 
pc
dd <- summary(pc) 


################################################################################
### Assuming you have PCA results stored in an object, e.g., pca_result
pc_df<-data.frame(fish_id=rownames(as.data.frame(pc$x)),pc$x)
write.csv(pc_df, file = "data/pca_scores.csv",row.names = FALSE)  # Save scores
#write.csv(pc$rotation, file = "pca_loadings2.csv")  # Save loadings
#write.csv(dd$importance, file = "pca_summary.csv")  # Save explained variance

################################################################################

#calculating mean pc scores and standard error for graphing later

PC1 <- tapply(pc$x[,1], info$fish_id, mean)
PC1

PC2 <- tapply(pc$x[,2], info$fish_id,mean)
PC2

PC3<- tapply(pc$x[,3], info$fish_id,mean)
PC3

plot(PC1,PC2)

plot(PC2,PC3)

combined_PC <- rbind(PC1, PC2, PC3)

combined_PC_t <- t(combined_PC)


library(ggplot2)

ggplot(combined_PC_t, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  labs(title = "",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

ggplot(combined_PC_t, aes(x = PC2, y = PC3)) +
  geom_point(size = 3) +
  labs(title = "",
       x = "PC2",
       y = "PC3") +
  theme_minimal()


st.err <- function(x) {+sd(x)/sqrt(length(x))}

PC1SE <- tapply(pc$x[,1], info$fish_id,st.err)

PC2SE <- tapply(pc$x[,2], info$fish_id,st.err)

PC3SE <- tapply(pc$x[,3], info$fish_id,st.err)
#site.pull

#mshape = Estimate the mean shape for a set of aligned specimens

Mean_morph <- mshape(gpa$coords)

GP <- gridPar(
  pt.bg = "red",
  pt.size = 1,
  link.col = "gray",
  link.lwd = 10,
  link.lty = 10,
  out.col = "gray",
  out.cex = 0.9,
  tar.pt.bg = "red",
  tar.pt.size = 1,
  tar.link.col = "pink",
  tar.link.lwd = 3,
  tar.link.lty = 1,
  tar.out.col = "black",
  tar.out.cex = 0.9,
  n.col.cell = 30,
  grid.col = "grey80",
  grid.lwd = 1,
  grid.lty = 1,
  txt.adj = NULL,
  txt.pos = 3,
  txt.cex = 0.8,
  txt.col = "black"
)


GP1 <- gridPar( #VECTOR
  pt.bg = "red",
  pt.size = 0.7,
  link.col = "gray",
  link.lwd = 10,
  link.lty = 10,
  out.col = "gray",
  out.cex = 0.9,
  tar.pt.bg = "red",
  tar.pt.size = 1,
  tar.link.col = "pink",
  tar.link.lwd = 3,
  tar.link.lty = 1,
  tar.out.col = "black",
  tar.out.cex = 0.9,
  n.col.cell = 30,
  grid.col = "grey80",
  grid.lwd = 1,
  grid.lty = 1,
  txt.adj = NULL,
  txt.pos = 3,
  txt.cex = 0.8,
  txt.col = "black"
)

PC1pred <- shape.predictor(gpa$coords, pc$x[,1], Intercept = FALSE, pred1=-0.04, pred2=0.04)

plotRefToTarget(Mean_morph, PC1pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC1pred$pred1, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC1pred$pred1, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)

plotRefToTarget(Mean_morph, PC1pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC1pred$pred2, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC1pred$pred2, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)


PC2pred <- shape.predictor(gpa$coords, pc$x[,2], Intercept = FALSE, pred1=-0.017, pred2=0.017)

plotRefToTarget(Mean_morph, PC2pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC2pred$pred1, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC2pred$pred1, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)

plotRefToTarget(Mean_morph, PC2pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC2pred$pred2, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC2pred$pred2, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)

PC3pred <- shape.predictor(gpa$coords, pc$x[,3], Intercept = FALSE, pred1=-0.011, pred2=0.011)

plotRefToTarget(Mean_morph, PC3pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC3pred$pred1, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC3pred$pred1, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)

plotRefToTarget(Mean_morph, PC3pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC3pred$pred2, mag = 6, method="points", gridPars=GP, label = TRUE)
plotRefToTarget(Mean_morph, PC3pred$pred2, mag = 4, method="vector", gridPars=GP1, label = TRUE, axes = TRUE, useRefPts = TRUE)


library(gridExtra)

# Create empty plots (plotRefToTarget doesn't return objects)
pdf("tables_figures/geo_plot1.pdf")
plotRefToTarget(Mean_morph, PC1pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()

pdf("tables_figures/geo_plot2.pdf")
plotRefToTarget(Mean_morph, PC1pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()

pdf("tables_figures/geo_plot3.pdf")
plotRefToTarget(Mean_morph, PC2pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()

pdf("tables_figures/geo_plot4.pdf")
plotRefToTarget(Mean_morph, PC2pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()

pdf("tables_figures/geo_plot5.pdf")
plotRefToTarget(Mean_morph, PC3pred$pred1, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()

pdf("tables_figures/geo_plot6.pdf")
plotRefToTarget(Mean_morph, PC3pred$pred2, mag = 6, method="TPS", gridPars=GP, label = TRUE)
dev.off()
