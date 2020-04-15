#########################################################################################
## 21
## Build the functional network
## Read the raster with patch ids 
## Compute patch's size and keep only those of size > 5ha
## These will become the nodes of the functional network
## Save a data frame with both ids (for patch and node), size, centroid's coordinates
#########################################################################################

build.fnet <- function(scn.path, scn.name, run, step, cluster, mrc){
  
  library(raster)
  library(FD)
  library(reshape2)
  library(igraph)
  library(tidyverse)
  
  if(!is.na(mrc))
    print(paste("Building functional network of", scn.name, "Replicate", run, "Step", step, "Clustering", cluster, "MRC", mrc))
  else
    print(paste("Building functional network of", scn.name, "Replicate", run, "Step", step, "Clustering", cluster))
  
  ###################### Load and read inputs and LANDIS raster outputs ######################
  
  ## Load the raster file (in .asc) with patch id of forest stands
  ## Load nodes of the functional network and links of the complete direct graph 
  if(is.na(mrc)){
    PATCH.ID <- raster(paste0("DataSp/ClusterID_", cluster, ".asc"))
    load(file=paste0("DataOut/NodesCDG_", cluster, ".rdata"))
    load(file=paste0("DataOut/LinksCDG_", cluster, ".rdata"))
  }
  else{
    PATCH.ID <- raster(paste0("DataSp/ClusterID_", cluster, "_mrc", mrc, ".asc"))
    load(file=paste0("DataOut/NodesCDG_", cluster, "_mrc", mrc, ".rdata"))
    load(file=paste0("DataOut/LinksCDG_", cluster, "_mrc", mrc, ".rdata"))
  }
  
  
  ## Dispersal capacity of tree species is Effective dispersal capacity taken from LANDIS
  spp.dispersal.dist <- read.table("DataIn/DispersalDist_19SppCdQ.txt", header=T) 
  
  ## Maturity age of tree species
  sex.maturity <- read.table("DataIn/SexMaturity_19sppCdQ.txt", header=T) 
  
  ## Species aggregated by functional group
  load(file="DataOut/FuncGroups.rdata")
  
  ## Load dissimilarity matrix
  load(file="DataOut/GowerDistanceGeneral.rdata")
  
  ## Build a RasterStack with the rasters accounting for the above ground biomass at year 0,
  ## build also the NSPP raster accounting for the number of species per cell (use this as a MASK)
  species <- c("abiebals", "acerrubr", "acersacc", "betualle", "betupapy", "betupopu",
               "fagugran", "fraxamer", "larilari", "piceabie", "piceglau", "picemari",
               "picerube", "pinuresi", "pinustro", "popugran", "poputrem", "thujocci", "tsugcana")
  AGB <- PATCH.ID
  names(AGB) <- "patch.id"
  
  for(spp in species){
    age.maturity <- sex.maturity$age[sex.maturity$species==spp]
    AGB.SPP <- raster( paste0(scn.path, sub("\\_.*", "", scn.name),"_rep", run,"/",scn.name, "/agbiomass/", spp, "/AGBiomass", step, ".img")  )
    AGE.SPP <- raster(paste0(scn.path, sub("\\_.*", "", scn.name),"_rep", run,"/",scn.name,"/MaxAge/", spp, "-", step, ".img"))
    AGB.SPP[] <- AGB.SPP[] * (AGE.SPP[]>=age.maturity) 
    names(AGB.SPP) <- spp
    extent(AGB.SPP) <- extent(PATCH.ID)
    AGB <- stack(AGB, AGB.SPP)
  }
  rm(AGB.SPP); rm(AGE.SPP); gc()
  
  
  ###################### Nodes's data frame informing above ground biomass, #######################
  ###################### functional group and 2 functional diversity indexes ######################
  
  ## Build a data frame with Above Ground Biomass per spp, and total AGB per ha
  spp.abund.node <- data.frame(patch.id=AGB$patch.id[],
                               abiebals=AGB$abiebals[], acerrubr=AGB$acerrubr[], acersacc=AGB$acersacc[],
                               betualle=AGB$betualle[], betupapy=AGB$betupapy[], betupopu=AGB$betupopu[],
                               fagugran=AGB$fagugran[], fraxamer=AGB$fraxamer[], larilari=AGB$larilari[],
                               piceabie=AGB$piceabie[], piceglau=AGB$piceglau[], picemari=AGB$picemari[],
                               picerube=AGB$picerube[], pinuresi=AGB$pinuresi[], pinustro=AGB$pinustro[],
                               popugran=AGB$popugran[], poputrem=AGB$poputrem[], thujocci=AGB$thujocci[],
                               tsugcana=AGB$tsugcana[])
  spp.abund.node <- left_join(spp.abund.node, nodesCDG)
  # 337.117 forest cells > 333.695 ha forest cells in the network
  spp.abund.node <- filter(spp.abund.node, !is.na(node.id)) %>% group_by(node.id, patch.id, size, xm, ym) %>%
    summarise(abiebals=sum(abiebals), acerrubr=sum(acerrubr), acersacc=sum(acersacc),
              betualle=sum(betualle), betupapy=sum(betupapy), betupopu=sum(betupopu),
              fagugran=sum(fagugran), fraxamer=sum(fraxamer), larilari=sum(larilari),
              piceabie=sum(piceabie), piceglau=sum(piceglau), picemari=sum(picemari),
              picerube=sum(picerube), pinuresi=sum(pinuresi), pinustro=sum(pinustro),
              popugran=sum(popugran), poputrem=sum(poputrem), thujocci=sum(thujocci),
              tsugcana=sum(tsugcana))
  
  ## Compute node's FUNCTIONAL DISPERSION with lingoes correction (Laliberte & Legendre 2010)
  ## Little tip: To avoid "fdisp" returning an error message when a patch (or community) 
  ## doesn't have any presence  --> At least one community has zero-sum abundances (no species).
  ## Assign 1 to one species (e.g. abiebals), so the function won't return any error and fdisp will be 0
  ## FDisp doesnt work  when a species has zero total abundance across all communities. 
  ## To avoid this, the following procedure is implemented below:  
  ## 1. look for a node with one single spp; add value 1 to that spp with zero abundance
  ## 2. compute fdisp for all the nodes 
  ## 3. change back to 0 the fdisp value for that node 
  ## 4. if there's no node with one single spp, then add the dummy value 1 to that missing spp, and  recompute fdisp for that node without the dummy value. 
  
  ################################################ TIP ############################################################
  ## Count number of spp present in each node
  spp.abund.node$nspp <- apply(spp.abund.node[6:24], 1, count.spp)
  
  ## Find spp without abundance for any spp
  spp.missing <- which(apply(spp.abund.node, 2, sum)==0)
  
  ## For communities without any spp, assign 1 to to one spp (e.g. abiebals)
  ## fdisp will be 0 anyways
  spp.abund.node$abiebals[spp.abund.node$nspp==0] <- 1
  
  ## For species not present in any node
  if(length(spp.missing)>0){
    if(any(spp.abund.node$nspp==1))
      node.id <- sample(spp.abund.node$node.id[spp.abund.node$nspp==1], 1, replace=F)
    else 
      node.id <- sample(spp.abund.node$node.id, 2, replace=F)
    spp.abund.node[spp.abund.node$node.id %in% node.id, spp.missing] <- 1  
  }
  
  ## Calculate functional dispersion with the modified matrix  
  fd <-  fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]) )
  
  ## In case there is a spp missing and no node with only one spp:
  if(length(spp.missing)>0 & !any(spp.abund.node$nspp==1)){
    ## compute true fd for node.id1
    spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 0
    fd.node.id1 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
    fd.node.id1 <- fd.node.id1$FDis[node.id[1]]
    # now for node.id2
    spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 1
    
    spp.abund.node[spp.abund.node$node.id == node.id[2], spp.missing] <- 0
    fd.node.id2 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
    fd.node.id2 <- fd.node.id2$FDis[node.id[2]]
    
    fd$FDis[node.id] <- c(fd.node.id1, fd.node.id2)
  }
  
  ## Back to original spp.abund.node
  spp.abund.node$abiebals[spp.abund.node$nspp==0] <- 0
  spp.abund.node[spp.abund.node$node.id %in% node.id, spp.missing] <- 0  
  
  ## Remove last colum
  spp.abund.node <- select(spp.abund.node, -nspp)
  
  ##############################################################################################################
  
  ## Compute node's functional diversity as Shannon index of functional groups
  aux <- melt(as.data.frame(spp.abund.node[,-(2:5)]), id=c("node.id"))
  names(aux)[2:3] <- c("species", "abund")
  aux <- left_join(aux, select(func.group, species, FG)) %>%
    group_by(node.id, FG) %>% summarize(abund=sum(abund))
  fg.abund.node <- acast(aux, node.id~FG)
  fg.abund.node <- data.frame(node.id=as.numeric(row.names(fg.abund.node)), fg.abund.node, 
                              fdiv=exp(vegan::diversity(fg.abund.node, "shannon")), nFG=apply(fg.abund.node>0,1,sum))
  names(fg.abund.node)[2:7] <- paste0("FG",1:6)
  
  ## Inform both diversity indexes in the main data frame 
  func.nodes <- as.data.frame(left_join(spp.abund.node, fg.abund.node) )
  func.nodes$fdisp <- fd$FDis
  
  
  ################################## Build the functional network ##################################
  
  ## 1. Build the 'spp.dist.node' data frame with node.id and 
  ## the dispersal capacity of the species present in each node
  ## Then, find the maximal dispersal distance per node, 'node.maxDD'
  spp.dist.node <- matrix(spp.dispersal.dist$d, nrow=nrow(func.nodes), 
                          ncol=nrow(spp.dispersal.dist), byrow=T)*(spp.abund.node[-(1:5)]>0) 
  spp.dist.node <- data.frame(node.id=func.nodes$node.id, spp.dist.node)
  node.maxDD <- data.frame(node.id=func.nodes$node.id, m=apply(spp.dist.node[-1], 1, max))
  
  ## 2. Keep the pair of nodes that the distance between them is lower than the maximal 
  ## dispersal distance of the source node. Joint the functional diversity index of the source node,
  ## and weight it by the number of functional groups represented in the source node (FG richeness).
  ## I do so to avoid situations in which less number of species, and in small quantity can
  ## travel to the sink node, but still have higher func.diversity (because of how the Shannon index is defined)
  ## e.g. from (70,10,10,10) to (10,10,0,10)
  ## Also join the functional dispersion index
  func.links <- left_join(linksCDG, node.maxDD, by=c("node1" = "node.id")) %>%  filter(dist<=m)  %>%
    left_join(select(func.nodes, node.id, fdiv, nFG, fdisp), by=c("node1" = "node.id")) %>% select(-m)
  func.links <- func.links[order(func.links$node1, func.links$node2),]  
  func.links$fdiv.node.w <- func.links$fdiv * pmax(log(func.links$nFG),0.5)
  names(func.links)[5:7] <- paste0(names(func.links)[5:7],".node")
  rm(node.maxDD); gc()
  
  ## 3. Add patch1 and patch2 corresponding ids
  func.links <- left_join(func.links, select(func.nodes, node.id, patch.id),  by=c("node1" = "node.id")) %>%
    left_join(select(func.nodes, node.id, patch.id),  by=c("node2" = "node.id"))
  names(func.links)[(ncol(func.links)-1):(ncol(func.links))] <- c("patch1", "patch2")
  
  ## 4. Build an auxiliar data.frame with the dispersal capacity of species in the source node of each link
  ## Then, mask species abundances of the source node that has lower dispersal capacity than the distance
  ## between it and the traget node
  spp.dist.link <- select(func.links, link.id, node1) %>%  left_join(spp.dist.node,  by=c("node1" = "node.id"))
  spp.abund.link <- select(func.links, link.id, node1) %>% left_join(spp.abund.node[,-(2:5)], by=c("node1" = "node.id"))
  spp.abund.link[-(1:2)] <- spp.abund.link[-(1:2)]  * (func.links$dist <= spp.dist.link[-(1:2)])
  
  ## 5. Compute link's weights following 3 criteria:
  ##   a. percentage of functional diversity 
  ##   b. percentage of weighted functional diversity 
  ##   c. percentage of functional dispersion
  ## that can travel from node source to node sink
  
  ## First, for each source node aggregate the species abundance in functional group abundances.
  ## For that, make explicit the correspondence between species and func.groups.
  ## Then join FG abundances per link to the "func.links" df
  aux <- melt(spp.abund.link, id=c("link.id", "node1")) %>% 
    left_join(select(func.group, species, FG), by=c("variable"="species"))
  aux <- group_by(aux, link.id, FG) %>% summarise(abund=sum(value))
  fg.abund.link <- data.frame(acast(aux, link.id~FG))
  fg.abund.link$link.id <- as.numeric(row.names(fg.abund.link))    
  func.links <- left_join(func.links, fg.abund.link)
  
  ## Second, compute the functional diversity (fdiv.link) of each link and fdiversity weighted
  func.links$fdiv.link <- exp(vegan::diversity(select(func.links, X1,X2,X3,X4,X5,X6), "shannon"))
  func.links$nFG.link  <- apply(select(func.links, X1,X2,X3,X4,X5,X6)>0, 1, sum)
  func.links$fdiv.link.w <- func.links$fdiv.link * pmax(log(func.links$nFG.link),0.5)
  func.links$pctg.fdiv <- pmin(func.links$fdiv.link/func.links$fdiv.node,1)
  func.links$pctg.fdiv.w <- pmin(func.links$fdiv.link.w/func.links$fdiv.node.w,1)
  
  ## Third, functional dispersion of links using the global dissimilartiy matrix
  ## TIP for "At least one community has zero-sum abundances (no species)."
  ## or "At least one species does not occur in any community (zero total abundance across all communities)."
  ################################################ TIP ############################################################
  ## Count number of spp present in each node
  spp.abund.node$nspp <- apply(spp.abund.node[6:24], 1, count.spp)
  
  ## Find spp without abundance for any spp
  spp.missing <- which(apply(spp.abund.node, 2, sum)==0)
  
  ## For communities without any spp, assign 1 to to one spp (e.g. abiebals)
  ## fdisp will be 0 anyways
  spp.abund.node$abiebals[spp.abund.node$nspp==0] <- 1
  
  ## For species not present in any node
  if(length(spp.missing)>0){
    if(any(spp.abund.node$nspp==1))
      node.id <- sample(spp.abund.node$node.id[spp.abund.node$nspp==1], 1, replace=F)
    else 
      node.id <- sample(spp.abund.node$node.id, 2, replace=F)
    spp.abund.node[spp.abund.node$node.id %in% node.id, spp.missing] <- 1  
  }
  
  ## Calculate functional dispersion with the modified matrix  
  aux <- rbind(spp.abund.link[,-c(1:2)], spp.abund.node[,-c(1:5,ncol(spp.abund.node))])
  fd <- fdisp(distrait, as.matrix(aux))
  fd <- fd[[1]]
  
  ## In case there is a spp missing and no node with only one spp:
  if(length(spp.missing)>0 & !any(spp.abund.node$nspp==1)){
    ## compute true fd for node.id1
    spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 0 
    fd.node.id1 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
    fd.node.id1 <- fd.node.id1$FDis[node.id[1]]
    # now for node.id2
    spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 1
    
    spp.abund.node[spp.abund.node$node.id == node.id[2], spp.missing] <- 0
    fd.node.id2 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
    fd.node.id2 <- fd.node.id2$FDis[node.id[2]]
    # fd$FDis[node.id] <- c(fd.node.id1, fd.node.id2)  #Marco 08.04.2020: I changed this line with the below one
    fd[node.id] <- c(fd.node.id1, fd.node.id2)
  }
  
  # ## Calculate functional dispersion with the modified matrix  
  # fd1 <-  fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]) )
  # 
  # ## In case there is a spp missing and no node with only one spp:
  # if(length(spp.missing)>0 & !any(spp.abund.node$nspp==1)){
  #   ## compute true fd for node.id1
  #   spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 0
  #   fd.node.id1 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
  #   fd.node.id1 <- fd.node.id1$FDis[node.id[1]]
  #   # now for node.id2
  #   spp.abund.node[spp.abund.node$node.id == node.id[1], spp.missing] <- 1
  #   
  #   spp.abund.node[spp.abund.node$node.id == node.id[2], spp.missing] <- 0
  #   fd.node.id2 <- fdisp(distrait, as.matrix(spp.abund.node[,-c(1:5,ncol(spp.abund.node))]))
  #   fd.node.id2 <- fd.node.id2$FDis[node.id[2]]
  #   
  #   fd$FDis[node.id] <- c(fd.node.id1, fd.node.id2)
  # }  
  
  
  
  ## Back to true spp.abund.node
  spp.abund.node$abiebals[spp.abund.node$nspp==0] <- 0
  spp.abund.node[spp.abund.node$node.id %in% node.id, spp.missing] <- 0  
  
  ## Remove last colum
  spp.abund.node <- select(spp.abund.node, -nspp)
  
  ##############################################################################################################
  
  ## Assign fd to links
  func.links$fdisp.link <- fd[1:nrow(spp.abund.link)]
  func.links$pctg.fdisp <- ifelse(func.links$fdisp.node==0, 0,
                                  pmin(func.links$fdisp.link/func.links$fdisp.node, 1))
  
  ## 6. Define a graph object using the two columns in func.links that define 
  ## the source - sink nodes of each link (called an edge list in igraph).
  ## It returns a 'igraph'
  ## By doing func.net[[i]], returns a list with the node.id of i's neighbours
  ## By doing func.net[i], returns a vector of length num.nodes with 0 and 1, masking i's neighbours
  func.links.fdisp <- func.links[func.links$pctg.fdisp>0,]
  func.net.fdisp <- graph(t(func.links.fdisp[,2:3]), directed=T)
  func.net.fdiv <- graph(t(func.links[,2:3]), directed=T)
  
  ## 7. Add vertex and edge attributes
  func.net.fdisp <- enrich.net(func.net.fdisp, func.nodes, func.links.fdisp)
  func.net.fdiv <- enrich.net(func.net.fdiv, func.nodes, func.links)
  
  ## 8. Save nodes and links of the functional network as well the network itself
  if(is.na(mrc)){
    save(func.nodes, file=paste0("DataScn/NodesFNET_", scn.name, "_run", run, "_step", step, "_", cluster, ".rdata"))
    save(func.links, file=paste0("DataScn/LinksFNET_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
    func.net <- func.net.fdisp
    save(func.net, file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
    func.net <- func.net.fdiv
    save(func.net, file=paste0("DataScn/FuncNetFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
  } else {
    save(func.nodes, file=paste0("DataScn/NodesFNET_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".rdata"))
    save(func.links, file=paste0("DataScn/LinksFNET_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc,".Rdata"))
    func.net <- func.net.fdisp
    save(func.net, file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster,"_mrc", mrc, ".Rdata"))
    func.net <- func.net.fdiv
    save(func.net, file=paste0("DataScn/FuncNetFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc,".Rdata"))
  }
  
  ## 9. Write tables of nodes and links to be used by CONEFOR, as many as 
  ## node.attributes - edge-attributes we aim to analyze
  if(is.na(mrc)){
    write.table(select(func.nodes, node.id, fdiv), paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.nodes, node.id, fdisp), paste0("DataScn/NodeFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdiv), paste0("DataScn/LinkFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdiv.w), paste0("DataScn/LinkFDIVW_", scn.name, "_run", run, "_step", step, "_", cluster, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdisp), paste0("DataScn/LinkFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".txt"), col.names=F, row.names=F, sep="\t")
  } else{
    write.table(select(func.nodes, node.id, fdiv), paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.nodes, node.id, fdisp), paste0("DataScn/NodeFDISP_", scn.name, "_run", run, "_step", step, "_", cluster,"_mrc", mrc, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdiv), paste0("DataScn/LinkFDIV_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdiv.w), paste0("DataScn/LinkFDIVW_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc,".txt"), col.names=F, row.names=F, sep="\t")
    write.table(select(func.links, node1, node2, pctg.fdisp), paste0("DataScn/LinkFDISP_", scn.name, "_run", run, "_step", step, "_", cluster,"_mrc", mrc, ".txt"), col.names=F, row.names=F, sep="\t")
  }
  
}


enrich.net <- function(func.net, func.nodes, func.links){
  # Node id (correlative id for the reduced network)
  V(func.net)$node.id <- func.nodes$node.id
  # Patch id (from the full landscape)
  V(func.net)$patch.id <- func.nodes$patch.id
  # Node area
  V(func.net)$size <- func.nodes$size
  # X and Y coordinates of the node's centroid
  V(func.net)$xm <- func.nodes$xm
  V(func.net)$ym <- func.nodes$ym
  # Node functional diversity, number of FG and functional dispersion 
  V(func.net)$fdiv <- func.nodes$fdiv
  V(func.net)$nFG <- func.nodes$nFG
  V(func.net)$fdisp <- func.nodes$fdisp
  # Node id of the source and sink nodes
  E(func.net)$node1 <- func.links$node1
  E(func.net)$node2 <- func.links$node2
  # Patch id of source and sink nodes
  E(func.net)$patch1 <- func.links$patch1
  E(func.net)$patch2 <- func.links$patch2
  # Edge length (Euclidian distance between nodes)
  E(func.net)$dist <- func.links$dist
  # Weight a. : percentage of functional diversity
  E(func.net)$pctg.fdiv <- func.links$pctg.fdiv
  E(func.net)$pctg.fdiv.w <- func.links$pctg.fdiv.w
  E(func.net)$pctg.fdisp <- func.links$pctg.fdisp
  return(func.net)
}


count.spp <- function(x){
  return(sum(x>0, na.rm = T))
}
