#########################################################################################
## 18
## Compute the dissimilarity matrix based on Gower's distance for those species in the landscape. 
## The  function 'gowdis' returns an object of class 'dist', i.e. a matrix 
## (nspp-1)x(nspp-1) whith NULL upper diagonal. 
## We will compute the generalization of the Gower's distance
#########################################################################################

rm(list=ls())
library(FD)
library(reshape2)
# setwd("C:/WORK/REDEFINE")
setwd("C:/Users/marco/OneDrive - UQAM/Documents/REDEFINE/Network/FunctionalNetwork/")

######################### BEFORE ADOPTING LANDIS VALUES, OCT2019 #########################
## Load functional traits of the species-
func.traits <- read.table("DataIn/FuncTraits_19SppCdQ.txt", header=T)
row.names(func.traits) <- func.traits$species

## Classifity traits by type
# Quantitative traits: tolerance related traits, wood density, max height, 
# seed mass, leaf mass area, Nmass
quant.trait <- func.traits[,3:10]
# Dichotomous traits: Association mycorhizienne (AM and ECM variables) 
dicho.trait <- data.frame(ifelse(func.traits[,11:12]=="N",0,1))
# Nominal traits: phylogenetic division and mode of reproduction a
nominal.trait <- data.frame(AM=func.traits[,13], reprod=as.factor(func.traits[,14]))
# Binary trait: mode of dispersion
binary.trait <- prep.binary(func.traits[,15:21], col.blocks = 7, label = "dispersal")
# Create a list of k data frames
ktab <- ktab.list.df(list(quant.trait, dicho.trait, nominal.trait, binary.trait))

## Now compute the dissimilarity Gowers distance matrix
distrait <- dist.ktab(ktab, c("Q", "D", "N", "B"), c("scaledBYrange"))
save(distrait, file="DataOut/GowerDistanceGeneral.rdata")



######################### WITH LANDIS VALUES, OCT2019 #########################
## Load functional traits of the species-
func.traits <- read.table("DataIn/FuncTraitsLANDIS_19SppCdQ_v2.txt", header=T)
row.names(func.traits) <- func.traits$species
func.traits
## Classifity traits by type
# Quantitative traits: tolerance related traits, wood density, max height, seed mass
quant.trait <- func.traits[,3:8]
# Nominal traits: phylogenetic division (G is conifers, A is deciduous)
# and main vector of seed dispersal: W is wind and A is animal
      # and mode of reproduction (1 is mostly vegetative, 2 is mostly by seed, 3 is only by seed)
      # nominal.trait <- data.frame(AG=func.traits[,9], reprod=as.factor(func.traits[,10]))
nominal.trait <- data.frame(AG=func.traits[,9], Disp=func.traits[,10])
# Create a list of k data frames
ktab <- ktab.list.df(list(quant.trait, nominal.trait))

## Now compute the dissimilarity Gowers distance matrix
distrait <- dist.ktab(ktab, c("Q", "N"), c("scaledBYrange"))
save(distrait, file="DataOut/GowerDistanceGeneral.rdata")
