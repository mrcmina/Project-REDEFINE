################################################################################
## 19
## Clustering tree species into functional groups using their functional traits
## Th e script was originally authored by Alain Paquette and then adapted by Kyle. T. Martins
## Current version is improved and customized to CdQ landscape
################################################################################

rm(list=ls())
library(FD)
library(cluster)
library(clValid)
library(vegan)
# setwd("C:/WORK/REDEFINE")
setwd("C:/Users/marco/OneDrive - UQAM/Documents/REDEFINE/Network/FunctionalNetwork/")

## Load species fuctional traits and the dissimilarty matrix
func.traits <- read.table("DataIn/FuncTraitsLANDIS_19SppCdQ_v2.txt", header=T)
load(file="DataOut/GowerDistanceGeneral.rdata")

## Conduct an agglomerative clustering on the Gower dissimilarity matrix of the 19 species
clust <- agnes(distrait, diss=TRUE, method = "ward")

## Research of the best number of groups (clusters)
msil <- NULL; mdunn <- NULL;   mconn <- NULL
nb <- 15
for (i in 2:nb) {
  classifcuts <- cutree(clust, k=i)
  silhouettes <- silhouette(classifcuts,  distrait)
  msil[i] <- mean(silhouettes[,3])
  mdunn[i] <- dunn(distrait,classifcuts)
  mconn[i] <- connectivity(distrait,classifcuts,neighbSize=10)
}

## Plot helpful graphs to identify the best number of groups
## Connectivity should be minimized, Dunn should be maximized as should silhouette
## useful link for interpretation of indices : 
## http://www.sthda.com/english/wiki/print.php?id=241
pdf("DataOut/FuncTraitsClusterIndexesLANDIS.pdf", width=8, height=4)
op <- par(mfrow=c(1,3))
plot(1:nb, msil,type="l", col="red", lwd=2, ylab="Average silhouette width", xlab="Number of groups")
for (k in 1:nb) 
  abline(v=k, col='orange', lty = "dotted")
grid(col = "orange", lty = "dotted", lwd = 1)
plot(1:nb, mdunn,type="l", col="blue", lwd=2, ylab="Dunn index", xlab="Number of groups")
for (k in 1:nb) 
  abline(v=k,col='orange',lty = "dotted")
grid(col = "orange", lty = "dotted", lwd = 1)
plot(1:nb, mconn,type="l", col="green", lwd=2, ylab="Connectivity index", xlab="Number of groups")
for (k in 1:nb) 
  abline(v=k,col='orange',lty = "dotted")
grid(col = "orange", lty = "dotted", lwd = 1)
par(op)
dev.off()
  
## Between 6 and 8 groupes would be the most valid; going to go with 6
pdf("DataOut/FuncTraitsClusterTreeLANDIS_6groups.pdf", width=8, height=4)
plot(clust, cex=0.5, which.plots=2, type="triangle")
# rect.hclust(clust, k=3, border="blue")
rect.hclust(clust, k=6, border="forestgreen") 
# rect.hclust(clust, k=7, border="red") 
dev.off()


## Plot the dendrogram of the 6 groups
# library(ggdendro); library(dendextend)
# plot(clust, cex=0.5, which.plots=2)
# rect.hclust(clust, k=6, border="forestgreen") 
# 
# clust_dendro <- as.dendrogram(clust)
# clust_dendro  %>% rect.dendrogram(k=6) %>% plot()



## Merge the functional groups with the functional group table
k <- cutree(clust, k=6)
func.group <- data.frame(func.traits, FG=k)
save(func.group, file="DataOut/FuncGroups.rdata")



    # ## Create a biplot to visualze and associate the traits with each functional group
    # NMDS1 <- metaMDS(distrait)
    # plot(NMDS1, type="t")
    # ordihull(NMDS1, groups=k, draw="polygon", col="grey90", label=F)
    # fit <- envfit(NMDS1, func.traits, perm=9999)
    # plot(fit, cex=1)
    # cent <- fit$factors$centroids
    # text(cent, labels=rownames(cent), col="red", cex=1.5)

   

