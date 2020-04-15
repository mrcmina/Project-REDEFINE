#########################################################################################
## 22
## Analyze a functional network, compute the following properties
## 1. Functional diversity
## 2. Functional redundancy
## 3. Modularity
## 4. Centrality
## 5. Connectivity
#########################################################################################

analyze.fnet <- function(scn.name, run, scn.step, cluster, mrc, fd="fdisp", node.importance=F, MM.PC=T){
  
  library(FD)
  library(moments)
  library(dplyr)
  library(igraph)
  
  if(is.na(mrc))
    print(paste("Analyzing functional network of", scn.name, "Step", scn.step, "Clustering", cluster))
  else
    print(paste("Analyzing functional network of", scn.name, "Step", scn.step, "Clustering", cluster, "MRC", mrc))
  
  ## Load functional network and its nodes
  if(is.na(mrc)){
    load(file=paste0("DataScn/NodesFNET_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".Rdata"))
    if(fd=="fdisp")
      load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".Rdata"))
    else
      load(file=paste0("DataScn/FuncNetFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".Rdata"))
  }
  else{
    load(file=paste0("DataScn/NodesFNET_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".Rdata"))
    if(fd=="fdisp")
      load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster,  "_mrc", mrc, ".Rdata"))
    else
      load(file=paste0("DataScn/FuncNetFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc,  ".Rdata"))  
  }
  
  ## FUNCTIONAL DIVERSITY
  print("01. Functional Diversity")
  fg.abund.net <- matrix(func.nodes$size, nrow=nrow(func.nodes), ncol=6, byrow=F) * 
                  (select(func.nodes, FG1, FG2, FG3, FG4, FG5, FG6)/1000) 
  fg.abund.net <- apply(fg.abund.net, 2, sum)/sum(func.nodes$size)
  net.fdiv <- exp(vegan::diversity(fg.abund.net))
  print(net.fdiv)
  
  
  ## FUNCTIONAL DISPERSION
  print("02. Functional Dispersion")
  load("DataOut/GowerDistanceGeneral.rdata")
  spp.abund.net <- matrix(func.nodes$size, nrow=nrow(func.nodes), ncol=19, byrow=F) * 
                   (select(func.nodes,  -patch.id, -size, -xm, -ym, -FG1, -FG2, -FG3, -FG4, -FG5, -FG6, -fdiv, -nFG, - fdisp)/1000) 
  
  ## Error-catcher for fdisp: add value 1 to nodes with 
  spp.abund.net$nspp <- apply(spp.abund.net, 1, count.spp) ## Count number of spp present in each node
  spp.missing <- which(apply(spp.abund.net, 2, sum)==0) ## Find spp without abundance for any spp
  
  ## For species not present in any node
  if(length(spp.missing)>0){
    if(any(spp.abund.net$nspp==1))
      node.id <- sample(spp.abund.net$node.id[spp.abund.net$nspp==1], 1, replace=F)
    else 
      node.id <- sample(spp.abund.net$node.id, 2, replace=F)
    spp.abund.net[spp.abund.net$node.id %in% node.id, spp.missing] <- 1  
  }
  
  spp.abund.net <- select(spp.abund.net, -node.id, -nspp)
  spp.abund.net <- apply(spp.abund.net, 2, sum)/sum(func.nodes$size) 
  spp.abund.net[spp.abund.net==0] <- 1e-12 #add a minimum value should there still be some zeros
  
  net.fdisp <- fdisp(distrait, t(spp.abund.net))[[1]]
  
  ## Back to original spp.abund.net, in case it was modified above
  spp.abund.net <- matrix(func.nodes$size, nrow=nrow(func.nodes), ncol=19, byrow=F) * 
    (select(func.nodes, -node.id, -patch.id, -size, -xm, -ym, -FG1, -FG2, -FG3, -FG4, -FG5, -FG6, -fdiv, -nFG, - fdisp)/1000) 
  spp.abund.net <- apply(spp.abund.net, 2, sum)/sum(func.nodes$size) 
  
  print(net.fdisp)
  
  
  ## FUNCTIONAL REDUNDANCY
  print("03.Functional Redundancy")
  RaoQ <- divc(as.data.frame(spp.abund.net), distrait, F)
  relabund.net <- spp.abund.net/sum(spp.abund.net)
  net.fredund<- 1-RaoQ/(1-sum(relabund.net^2))
  print(as.numeric(net.fredund))
  
  
  # ## MODULARITY
  print("04.Modularity")
  # Community detection methods:
  # cluster_walktrap, cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass 
    # 2. Because the cluster_edge_betweenness performs an algorithm that calculates first the edge betweenness 
    # of each edge on the graph, and then removes the edge with the highest edge betweenness score, 
    # to then recalculate edge betweenness of the edges and again removing the one with the highest score, etc.
    # It takes an eternity. So, dismiss.
    # 3. Error in cluster_fast_greedy
    # fast greedy community detection works for undirected graphs only, Unimplemented function call
    # 4. cluster_spinglass 
    # The input graph, can be directed but the direction of the edges is neglected.
  if(fd=="fdisp")
    clust.wt <- cluster_walktrap(func.net, weights=E(func.net)$pctg.fdisp, steps=5, merges=T, modularity=T, membership=T)
  else
    clust.wt <- cluster_walktrap(func.net, weights=E(func.net)$pctg.fdiv, steps=5, merges=T, modularity=T, membership=T)
  modularity.wt <- modularity(clust.wt)
  print(modularity.wt)   
      # modularity.wt <- numeric(9)
      # for(i in 2:10){
      #   print(i)
      #   clust.wt <- cluster_walktrap(func.net, weights=E(func.net)$pctg.fdisp, steps=i, merges=T, modularity=T, membership=T)
      #   modularity.wt[i] <- modularity(clust.wt)
      # }
      # plot(x=2:10, modularity.wt[-1])
      
  ## CONNECTIVITY
  print("05.Connectivity and 06.Centrality")
  # Create a directory to call conefor 
  if(!dir.exists(paste0("Conefor/", scn.name, "_", fd)))
    dir.create(paste0("Conefor/", scn.name, "_", fd))
  
  # Copy the node and link files in this directory 
  if(fd=="fdisp"){
    if(is.na(mrc)){
      file.copy(paste0("DataScn/NodeFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".txt"),
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster,".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T)  
    }
    else{
      file.copy(paste0("DataScn/NodeFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"),
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDISP_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T) 
    }
  }
  if(fd=="fdiv"){
    if(is.na(mrc)){
      file.copy(paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster,".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T)  
    }
    else{
      file.copy(paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T)
    }
  }
  if(fd=="fdivw"){
    if(is.na(mrc)){
      file.copy(paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDIVW_", scn.name, "_run", run, "_step", scn.step, "_", cluster, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T)  
    }
    else{
      file.copy(paste0("DataScn/NodeFDIV_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/NodeFD.txt"), overwrite=T)
      file.copy(paste0("DataScn/LinkFDIVW_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, ".txt"), 
                paste0("Conefor/", scn.name, "_", fd, "/LinkFD.txt"), overwrite=T)  
    
    }
  }
  # The conefor.exe too
  file.copy("Conefor/coneforWin64.exe", paste0("Conefor/", scn.name, "_", fd))
  
  # Change the working directory to this new folder
  if(MM.PC)
    setwd(paste0("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork/Conefor/", scn.name, "_", fd))  ## change!!!
  else
    setwd(paste0("C:/work/REDEFINE/Conefor/", scn.name, "_", fd))
  
  # Call conefor to calculate the global indexes and node importance if required
  conefor.cmd <- ifelse(node.importance,
                        paste("coneforWin64.exe -nodeFile", "NodeFD.txt", "-conFile", "LinkFD.txt", "-t prob notall -PC -BCPC"),
                        paste("coneforWin64.exe -nodeFile", "NodeFD.txt", "-conFile", "LinkFD.txt", "-t prob notall -PC onlyoverall") )
  print(conefor.cmd)
  shell(conefor.cmd)
  
  # Read the target results 
  conefor.result <- read.table("overall_indices.txt", header=F, sep="\t")
  names(conefor.result) <- c("index", "val")
  PCnum <- conefor.result$val[conefor.result$index=="PCnum"]
  EC_PC <- conefor.result$val[conefor.result$index=="EC(PC)"]
  if(fd=="fdisp")
    PC <- PCnum / (0.4*nrow(func.nodes))^2
  else
    PC <- PCnum / (6*nrow(func.nodes))^2
  if(node.importance){
    # Compute betweenness centrality
    aux <- read.table("node_importances.txt", header=T, sep="\t")
    aux <- aux[order(aux$Node),]
    centrality <- data.frame(btw=betweenness(func.net, V(func.net), directed=TRUE, weights=E(func.net)$pctg.fd.w,
                                                    nobigint=TRUE, normalized=F),
                             dPC=aux$dPC, dPCintra=aux$dPCintra,	dPCflux=aux$dPCflux, 
                             dPCconnector=aux$dPCconnector, dBC=aux$dBC_PC)
    fractions <- data.frame(component=c("intra", "flux", "connector"), 
                            teta=c(sum(aux$dPCintra)/sum(aux$dPC),
                                   sum(aux$dPCflux)/sum(aux$dPC),
                                   sum(aux$dPCconnector)/sum(aux$dPC)))
    dBCmoments <- data.frame(moment=c("mean", "variance", "median", "kurtosis", "skewness"),
                             val=c(mean(aux$dBC_PC), var(aux$dBC_PC), median(aux$dBC_PC),
                                   kurtosis(aux$dBC_PC), skewness(aux$dBC_PC)))
  }
  else{
    centrality <- data.frame(btw=NA, dPC=NA, dPCintra=NA, dPCflux=NA, dPCconnector=NA, dBC=NA)
    fractions <- data.frame(component=c("intra", "flux", "connector"), teta=c(NA,NA,NA))
    dBCmoments <- data.frame(moment=c("mean", "variance", "median", "kurtosis", "skewness"),
                             val=c(NA,NA,NA,NA,NA))
  }
  
  # Remove conefor.exe to save space and go back to the original working directory
  file.remove("coneforWin64.exe")
  if(MM.PC)
    setwd("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork")  ## change!!!
  else
    setwd("C:/work/REDEFINE")
  
  # Save the properties of the network
  net.prop <- list(net.fdiv=net.fdiv, net.fdisp=net.fdisp, net.fredund=as.numeric(net.fredund),
                   net.modularity=modularity.wt, dBCmoments=dBCmoments, fractions=fractions, 
                   centrality=centrality, PCnum=PCnum, EC=EC_PC, PC=PC)
  if(is.na(mrc))
    save(net.prop, file=paste0("DataScn/NetProp_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_", fd, ".Rdata"))
  else
    save(net.prop, file=paste0("DataScn/NetProp_", scn.name, "_run", run, "_step", scn.step, "_", cluster, "_mrc", mrc, "_", fd, ".Rdata"))
}

count.spp <- function(x){
  return(sum(x>0, na.rm = T))
}
