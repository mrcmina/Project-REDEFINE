#########################################################################################
## 17
## 27/09/2019
## Still a forest network need to be defined for each sub-region in CdQ
## Goal: networks with a reasonable number of nodes (<4000 ?) and none of them too big (<1000 ha?)
## The approaches and classifications tested up to now do not return satisfactory results
## Here, I first select the biggest patches to deal with
## Then I propose different forest type classification schemes to split them using both
## Above Ground Biomass (AGB) and a binary representation (i.e. presence) of AGB.
## Then I realize that a possible good approach is to deal patch by patch and use 
## coordinates x-y as a variable in the clustering process.
#########################################################################################

rm(list=ls())
library(reshape2)
library(gridExtra)
library(raster)
library(tidyverse)
library(landscapemetrics)
library(RStoolbox)
library(RANN)
setwd("C:/WORK/REDEFINE")


## Names of the species
species <- c("abiebals", "acerrubr", "acersacc", "betualle", "betupapy", "betupopu",
             "fagugran", "fraxamer", "larilari", "piceabie", "piceglau", "picemari",
             "picerube", "pinuresi", "pinustro", "popugran", "poputrem", "thujocci", "tsugcana")

## Scenario of reference
scn <- "LANDIS/CdQsims/Output/BAU_rcp45/rep1/agbiomass/"

## Read Initial communities raster file (with the right projection and extent) and the 
## shape file with the Management units, then rasterize it
FOREST <- raster("LANDIS/CdQsims/Input_Data/EcoregionsIC/initial_communities.tif")
if(0>1){
  MRC <- shapefile("C:/WORK/REDEFINE/DataIn/MRC_CdQ_MTM8/MRC_CdQ_MTM8.shp")
  rMRC <- rasterize(MRC, FOREST)
  writeRaster(rMRC, "DataSp/ManagementUnits.asc", format="ascii", NAflag=0, overwrite=T)
}
rMRC <- raster("DataSp/ManagementUnits.asc") 
## Look for NA in MRC that are not in FOREST and phagocite them, .
sum(is.na(rMRC[!is.na(FOREST[])]))
dta <- data.frame(mrc=rMRC[], frst=FOREST[], coordinates(rMRC))
nas <- filter(dta, is.na(mrc) & !is.na(frst))
size.neighs <- c(3^2, 5^2, 7^2, 9^2, 11^2)
step <- 1
while(nrow(nas)>0){
  neighs <- nn2(coordinates(rMRC), query=nas[,3:4], k=size.neighs[step], treetype="kd", searchtype="priority",  eps=0)
  for(i in 1:nrow(nas)){
    new.mrc <- round(median(dta[neighs$nn.idx[i,],1], na.rm=T))
    if(!is.na(new.mrc)){
      dta[neighs$nn.idx[i,1],1] <- new.mrc
      nas <- nas[-1,]
    }
  }
  step <- step +1
}
## Transfer new values to the raster
rMRC[] <- dta$mrc
plot(rMRC, col=rainbow(5))
table(rMRC[])
rm(nas); rm(size.neighs); rm(step); rm(new.mrc); rm(neighs); rm(dta); rm(i)
      # 1      2      3      4      5 
      # 118807 162753 123476 130369 191216


######################## MAKE ABOVE GROUND BIOMASS SPATIAL DATA READY ####################################
## Build a RasterStack with the rasters accounting for the above ground biomass at year 0,
## build also the NSPP raster accounting for the number of species per cell (use this as a MASK)
rm(AGB)
for(spp in species){
  AUX <- raster(paste0(scn, "/", spp, "/AGBiomass0.img"))
  names(AUX) <- spp
  if(exists("AGB"))
    AGB <- stack(AGB, AUX)
  else
    AGB <- AUX
}
rm(AUX); rm(spp);  gc()


## Build a data frame with Above Ground Biomass per spp, and total AGB per ha
agb.spp <- data.frame(abiebals=AGB$abiebals[], acerrubr=AGB$acerrubr[], acersacc=AGB$acersacc[],
                      betualle=AGB$betualle[], betupapy=AGB$betupapy[], betupopu=AGB$betupopu[],
                      fagugran=AGB$fagugran[], fraxamer=AGB$fraxamer[], larilari=AGB$larilari[],
                      piceabie=AGB$piceabie[], piceglau=AGB$piceglau[], picemari=AGB$picemari[],
                      picerube=AGB$picerube[], pinuresi=AGB$pinuresi[], pinustro=AGB$pinustro[],
                      popugran=AGB$popugran[], poputrem=AGB$poputrem[], thujocci=AGB$thujocci[],
                      tsugcana=AGB$tsugcana[])
agb.spp$tot <- apply(agb.spp,1,sum) # 1.649.532
agb.spp$mrc <- rMRC[]
agb.spp <- filter(agb.spp, tot>0)   # 337.091



########################### SPATIALLY DELIMITATE FOREST PATCHES ###########################
################ THEN USE A FOREST TYPE CLASSIFICATION TO SPIT THE BIGGEST ################
for(reg in 1:5){
  print(paste("Region", reg))
  CLASS <- FOREST
  CLASS[!is.na(CLASS[])] <-  1*(agb.spp$tot>0 & agb.spp$mrc==reg)
  CLASS[CLASS[]==0] <- NA
  PATCH.ID <- get_patches(CLASS, directions=8)[[1]]
  if(reg==1){
    patches <- data.frame(id=PATCH.ID[], coordinates(PATCH.ID))
    patches$reg <- NA
  }
  else{
    aux <- data.frame(id2=PATCH.ID[], coordinates(PATCH.ID)) 
    aux$id2 <- aux$id2 + max(patches$id, na.rm=T)
    aux <- left_join(patches, aux)
    subst <- which(!is.na(aux$id2) & is.na(aux$id))
    patches$id[subst] <- aux$id2[subst]
    patches$reg[subst] <- reg
  }
}
patches.area <- filter(patches, !is.na(id)) %>% group_by(id) %>% 
                summarize(area=length(id), xm=mean(x), ym=mean(y), reg=mean(reg))
nrow(patches.area)
#  2610

## The PATCH.ID raster
PATCH.ID[] <- patches$id
writeRaster(PATCH.ID, "DataSp/PatchID.asc", format="ascii", NAflag=-1, overwrite=T)  


## Number of patches with area <= small.th, and the percentage these represent
small.th <- 5
big.th <- 1500
nnode <- nrow(filter(patches.area, area<=small.th))
nnode; round(nnode/nrow(patches.area)*100,1)
# [1] 1624
# [1] 62.4
## Number of patches with area > small.th, and the percentage these represent
nnode <- nrow(patches.area) - nrow(filter(patches.area, area<=small.th))
nnode; round(nnode/nrow(patches.area)*100,1)
# [1] 986
# [1] 37.8
## Number of patches with area > big.th, and the percentage these represent
nnode <- nrow(filter(patches.area, area>big.th))
nnode; round(nnode/nrow(patches.area)*100,2)
# [1] 17
# [1] 0.65
## Number of patches with small.th < area <= big.th, and the percentage these represent
nnode <- nrow(filter(patches.area, area>small.th, area<=big.th))
nnode; round(nnode/nrow(patches.area)*100,1)
# [1] 969
# [1] 37.1


###################################### MANUALLY SPLIT THE BIGGEST ONES ######################################
## Before doing any clustering, manually split thse patch in a few sub patches by 
## artifically removing some bridging cells  (i.e. change their type to non-forest)
## Then apply the current clustering by using the coordinates - forest type approach.
## First, round coordinates in the "patches" data frame
biggest <- data.frame(indx=1:nrow(patches), round(coordinates(PATCH.ID)))
cell.indx <- numeric(0)
## Second, change to NA the target cells and reassing these values to the raster
# 01 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==396330 & biggest$y==5076263])
# 02 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==456430 & biggest$y==5136263])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==450030 & biggest$y==5133863])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==450030 & biggest$y==5133763])
# 03 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==416530 & biggest$y==5154663])
# 04 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==422230 & biggest$y==5149663])
# 05 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==393930 & biggest$y==5115863])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==394030 & biggest$y==5115863])
# 06 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==387930 & biggest$y==5112663])
# 07 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==384030 & biggest$y==5111463])
# 08 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==369730 & biggest$y==5084863])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==369830 & biggest$y==5084763])
# 09 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==395630 & biggest$y==5106363])
# 10 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==402230 & biggest$y==5108063])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==402830 & biggest$y==5106863])
# 11 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==434130 & biggest$y==5090863])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==437530 & biggest$y==5087263])
# 12 split 
cell.indx <- c(cell.indx, biggest$indx[biggest$x==450630 & biggest$y==5144963])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==448430 & biggest$y==5142663])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==447030 & biggest$y==5141263])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==445630 & biggest$y==5139763])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==443530 & biggest$y==5137563])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==442830 & biggest$y==5136863])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==440330 & biggest$y==5134463])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==438530 & biggest$y==5132363])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==436430 & biggest$y==5130263])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==434330 & biggest$y==5127463])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==431530 & biggest$y==5124463])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==430330 & biggest$y==5123163])
cell.indx <- c(cell.indx, biggest$indx[biggest$x==428630 & biggest$y==5121363])

## Third, import these changes to the CLASS raster before splitting patches by the spatial criteria
for(reg in 1:5){
  print(paste("Region", reg))
  CLASS <- FOREST
  class <- CLASS[]
  class[!is.na(class)] <-  1*(agb.spp$tot>0 & agb.spp$mrc==reg)
  class[class==0] <- NA
  class[cell.indx] <- NA
  CLASS[] <- class
  PATCH.ID <- get_patches(CLASS, directions=8)[[1]]
  if(reg==1){
    patches <- data.frame(id=PATCH.ID[], coordinates(PATCH.ID))
    patches$reg <- NA
    patches$reg[!is.na(patches$id)] <- reg
  }
  else{
    aux <- data.frame(id2=PATCH.ID[], coordinates(PATCH.ID)) 
    aux$id2 <- aux$id2 + max(patches$id, na.rm=T)
    aux <- left_join(patches, aux)
    subst <- which(!is.na(aux$id2) & is.na(aux$id))
    patches$id[subst] <- aux$id2[subst]
    patches$reg[subst] <- reg
  }
}
patches.area <- filter(patches, !is.na(id)) %>% group_by(id) %>% 
                summarize(area=length(id), xm=mean(x), ym=mean(y), reg=mean(reg))
nrow(patches.area)
sum(patches.area$area>1500)
# 2610  -> 2625 After splitting
# 17 > 26 After splitting



###################################### CLUSTERING OF BIGGEST PATCHES: COORDINATES  & FOREST TYPE ####################################
AGB.CONIF <- AGB$abiebals + AGB$piceabie + AGB$piceglau + AGB$picemari + AGB$picerube + AGB$larilari + AGB$thujocci + AGB$tsugcana + AGB$pinuresi + AGB$pinustro
AGB.DECID <- AGB$acerrubr + AGB$acersacc + AGB$fagugran + AGB$fraxamer + AGB$betualle + AGB$betupapy + AGB$betupopu + AGB$popugran + AGB$poputrem
CONIF <- FOREST; CONIF[] <- AGB.CONIF[]
DECID <- FOREST; DECID[] <- AGB.DECID[]

## Clustering of x-y coordinates of patches of size > 1500 by number of classes as patch.size/1000 ha
big.th <- 1500
last.id <- max(patches$id, na.rm=T)
patches <- left_join(patches, select(patches.area, id, area))  %>%  mutate(mask=ifelse(area>big.th, 1, NA))
idsbig <- sort(unique(patches$id[!is.na(patches$mask)]))
cells <- xyFromCell(FOREST, 1:ncell(FOREST))
COORD.X <- FOREST; COORD.X[] <- cells[,1] 
COORD.Y <- FOREST; COORD.Y[] <- cells[,2]
patches$clustidCoord <- NA
patches$clustidDominant <- NA
for(id.1big in idsbig){
  patches$mask.1big <- ifelse(patches$id==id.1big, 1, NA)
  MASK.1BIG <- FOREST; MASK.1BIG[] <- patches$mask.1big
  nClass <- max(2,round(sum(patches$mask.1big, na.rm=T)/1000))
  print(paste("ID:", id.1big, " Area:", patches.area$area[patches.area$id==id.1big], " NumClass:", nClass))
  set.seed(13)
  classCoord <- unsuperClass(brick(COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG), 
                             nClasses=nClass, nStarts=nClass*4, nIter=nClass*100, norm=T, clusterMap=F, algorithm = "Lloyd")
      # classAbund <- unsuperClass(brick(COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG,COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG, DECID*MASK.1BIG, CONIF*MASK.1BIG), 
      #                            nClasses=nClass, nStarts=nClass*4, nIter=nClass*100, norm=T, clusterMap=F, algorithm = "Lloyd")
      # classPresent <- unsuperClass(brick(COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG,COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG, (DECID>0)*MASK.1BIG, (CONIF>0)*MASK.1BIG), 
      #                              nClasses=nClass, nStarts=nClass*4, nIter=nClass*100, norm=T, clusterMap=F, algorithm = "Lloyd")
  classDominant <- unsuperClass(brick(COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG,COORD.X*MASK.1BIG, COORD.Y*MASK.1BIG, (DECID>=CONIF)*MASK.1BIG), 
                               nClasses=nClass, nStarts=nClass*4, nIter=nClass*100, norm=T, clusterMap=F, algorithm = "Lloyd")
  # plot(classDominant$map, col=rainbow(nClass))
  # writeRaster(classCoord$map, file=paste0("DataSp/ClusterCoord.PatchID", id.1big, ".nC", nClass, ".asc"), format="ascii", overwrite=T, NAflag=0)  
  #     # writeRaster(classAbund$map, file=paste0("DataSp/ClusterCoordAbund.PatchID", id.1big, ".nC", nClass, ".asc"), format="ascii", overwrite=T, NAflag=0)
  #     # writeRaster(classPresent$map, file=paste0("DataSp/ClusterCoordPresent.PatchID", id.1big, ".nC", nClass, ".asc"), format="ascii", overwrite=T, NAflag=0)
  # writeRaster(classDominant$map, file=paste0("DataSp/ClusterCoordDominant.PatchID", id.1big, ".nC", nClass, ".asc"), format="ascii", overwrite=T, NAflag=0)
  patches$clustidCoord[!is.na(patches$mask.1big) & patches$mask.1big==1] <- classCoord$map[patches$mask.1big==1]+last.id
  patches$clustidDominant[!is.na(patches$mask.1big) & patches$mask.1big ==1] <- classDominant$map[patches$mask.1big ==1]+last.id
  last.id <- last.id+nClass
}
## Transfer original id of small patches to new.id
patches$clustidCoord[is.na(patches$clustidCoord)] <- patches$id[is.na(patches$clustidCoord)]
patches$clustidDominant[is.na(patches$clustidDominant)] <- patches$id[is.na(patches$clustidDominant)]

## Remove inecessary columns in patches
patches <- select(patches, id, x, y, reg, clustidCoord, clustidDominant)

## Assing a cluster.id to the hidden cells
hidden <- patches[cell.indx,]
patches$xy <- paste0(patches$x, "-", patches$y)
for(i in 1:nrow(hidden)){
  a <- c(paste0(hidden$x[i]-100, "-", hidden$y[i]+100), 
         paste0(hidden$x[i], "-", hidden$y[i]+100),
         paste0(hidden$x[i]+100, "-", hidden$y[i]+100),
         paste0(hidden$x[i]-100, "-", hidden$y[i]),
         paste0(hidden$x[i]+100, "-", hidden$y[i]),
         paste0(hidden$x[i]-100, "-", hidden$y[i]-100),
         paste0(hidden$x[i], "-", hidden$y[i]-100),
         paste0(hidden$x[i]+100, "-", hidden$y[i]-100))
  b <- filter(patches, xy %in% a) %>% select(clustidCoord, clustidDominant)
  hidden$clustidCoord[i] <- round(median(b$clustidCoord, na.rm = T))
  hidden$clustidDominant[i] <- round(median(b$clustidDominant, na.rm=T))
}
patches$clustidCoord[cell.indx] <- hidden$clustidCoord
patches$clustidDominant[cell.indx] <- hidden$clustidDominant
patches <- select(patches, -xy)

## Checking  --> OK !!
    # > sum(!is.na(patches$clustidDominant))
    # [1] 337091
    # > sum(!is.na(patches$clustidCoord))
    # [1] 337091
    # > sum(!is.na(patches$id))
    # [1] 337061

## Write cluster IDs
PATCH.ID[] <- patches$clustidCoord
writeRaster(PATCH.ID, "DataSp/ClusterID_Coord.asc", format="ascii", NAflag=-1, overwrite=T)
PATCH.ID[] <- patches$clustidDominant
writeRaster(PATCH.ID, "DataSp/ClusterID_Dominant.asc", format="ascii", NAflag=-1, overwrite=T)

## Write cluster IDs by MRC
for(mrc in 1:5){
  aux <- patches$clustidCoord
  aux[!is.na(patches$reg) & patches$reg!=mrc] <- NA
  PATCH.ID[] <- aux
  writeRaster(PATCH.ID, paste0("DataSp/ClusterID_Coord_mrc", mrc, ".asc"), format="ascii", NAflag=-1, overwrite=T)
  aux <- patches$clustidDominant
  aux[!is.na(patches$reg) & patches$reg!=mrc] <- NA
  PATCH.ID[] <- aux
  writeRaster(PATCH.ID, paste0("DataSp/ClusterID_Dominant_mrc", mrc, ".asc"), format="ascii", NAflag=-1, overwrite=T)
  
}
