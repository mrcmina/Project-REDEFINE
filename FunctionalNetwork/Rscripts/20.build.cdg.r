#########################################################################################
## 20
## Build the complete direct graph for the CdQ landscape
## Read the raster with patch ids 
## Compute patch's size and keep only those of size > 5ha
## These will become the nodes of the functional network
## Save a data frame with both ids (for patch and node), size, centroid's coordinates
#########################################################################################

rm(list=ls())
library(raster)
library(FD)
library(reshape2)
library(tidyverse)
setwd("C:/Users/marco/OneDrive - UQAM/Documents/REDEFINE/Network/FunctionalNetwork/") # setwd("C:/WORK/REDEFINE")

## Mode of clustering
for(clustering in c("Dominant", "Coord")){
  for(mrc in c(NA, 1:5)){
    
    print(paste("Building Complete Directe Graph for MRC", mrc, "and clustering approach", clustering))
    ## Function to compute euclidean distance between patches from 
    ## their X - Y coordinates
    euclid.dist <- function(x, y) min(as.numeric(spDists(as.matrix(x), as.matrix(y))))
    
    ## Load the raster file (in .asc) with patch id of forest stands. Build a data frame 
    ## with patch id and cells' XY coordinates. Raster resolution is 1 ha
    ## Raster comes from 17.ftype.net.cluster.r script
    if(is.na(mrc))
      PATCH.ID <- raster(paste0("DataSp/ClusterID_", clustering, ".asc"))
    if(!is.na(mrc))
      PATCH.ID <- raster(paste0("DataSp/ClusterID_", clustering, "_mrc", mrc, ".asc"))
    patch.id <- data.frame(PATCH.ID[], coordinates(PATCH.ID))
    names(patch.id)[1] <- "id"
    
    ## Build a data frame with node.id, patch.id, area, and centroid coordinates
    small.th <- 5
    nodesCDG <- filter(patch.id, !is.na(id)) %>% group_by(id) %>% 
      summarise(size=length(id), xm=mean(x), ym=mean(y)) %>% filter(size > small.th) %>%
      mutate(patch.id=id) %>% select(-id)
    nodesCDG$node.id <- 1:nrow(nodesCDG)
    nodesCDG <- as.data.frame(nodesCDG[,c(5,4,1:3)])
    if(is.na(mrc))
    save(nodesCDG, file=paste0("DataOut/NodesCDG_", clustering, ".rdata"))
    if(!is.na(mrc))
      save(nodesCDG, file=paste0("DataOut/NodesCDG_", clustering, "_mrc", mrc,".rdata"))
    
    ## Compute the Euclidean distance between the nodes of the network, i.e.:
    ## For each pair of nodes save in a matrix: link.id, node1, node2, dist.from1.to2
    ## I want to calculate the distance from each target node to all the other nodes
    ## except for those nodes that have already been target nodes (to not repeat calculation)
    ## So, I will calculate n-1, n-2, n-3, .... distances. The sum of n integers is (n*(n + 1))/2
    ## I first work with patch.id as these are the ones informed in PATCH.ID raster
    nnode <- nrow(nodesCDG)
    ids <- nodesCDG$patch.id
    # List with X - Y coordinates of all cells within nodes
    coord <- vector("list", nnode)
    for(i in 1:nnode)
      coord[[i]] <- filter(patch.id, id==ids[i]) %>% select(-id)
    # Matrix to store the links, i.e. distance between pair of patches
    dist.patch <- matrix(nrow=((nnode-1)*nnode)/2, ncol=4)
    dist.patch[,1] <- 1:nrow(dist.patch)
    index.ini <- 1; index.end <- nnode-1
    for(i in 1:(nnode-1)){
      print(paste0("node: ", i,  "/", nnode, ", ", round(100*i/nnode,1), "%"))
      dist.patch[index.ini:index.end,2] <- ids[i]
      dist.patch[index.ini:index.end,3] <- ids[(i+1):nnode]
      dist.patch[index.ini:index.end,4] <- 
        sapply(coord[(i+1):nnode], euclid.dist, y = coord[[i]], simplify = "array")   
      index.ini <- index.end+1
      index.end <- index.ini + length((i+2):nnode) -1
    }
    
    ## Convert to a data frame and remplace patch.id by node.id
    dist.patch <- as.data.frame(dist.patch)
    names(dist.patch) <- c("link.id", "patch1", "patch2", "dist")
    linksCDG <- left_join(dist.patch, nodesCDG[,1:2], by=c("patch1"="patch.id")) %>%
      left_join(nodesCDG[,1:2], by=c("patch2"="patch.id")) 
    linksCDG <- linksCDG[, c(1,5,6,4)]
    names(linksCDG) <- c("link.id", "node1", "node2", "dist")
    
    ## Make the graph directed by copying pairs of source-target nodes, changing their role but assigning the same distance
    linksCDG <- rbind(linksCDG,
                      data.frame(link.id=linksCDG$link.id+max(linksCDG$link.id), 
                                 node1=linksCDG$node2, node2=linksCDG$node1, dist=linksCDG$dist))
    
    ## Save it
    if(is.na(mrc))
      save(linksCDG, file=paste0("DataOut/LinksCDG_", clustering, ".rdata"))
    if(!is.na(mrc))    
      save(linksCDG, file=paste0("DataOut/LinksCDG_", clustering, "_mrc", mrc, ".rdata"))
    
  }
}
