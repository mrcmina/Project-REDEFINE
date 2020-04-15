#########################################################################################
## 22
## Run the function that analyze the functional network and compute properties (fdiv, fred, Mod, Centr, Conn)
## It saves one .Rdata with network properties per scenario in folder DataScn 
#########################################################################################

rm(list=ls());gc()  #setwd("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork")
source("Rscripts/22.analyze.fnet.r")

## Analyze a set of scenarios every time step
mgmt.scenarios <-  c("NM", "BAU", "BAUDistHigh")

scenarios <- c(paste0(mgmt.scenarios,"_", "current"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp45"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp85"))


for(scn.name in scenarios){
  for(mrc in c(NA, 1:5)){ 
    for(run in 1:5){
      for(scn.step in seq(0,90,10))
        
        analyze.fnet(scn.name, run, scn.step, cluster="Dominant", mrc, fd="fdisp", node.importance=F, MM.PC=T)
    }
  }
}



## TEST Run the function that analyze the functional network, it saves a .Rdata with network properties
# scn.name <- "NM_current"
# scn.step <- 90
# cluster <- "Dominant"
# mrc <- NA
# run <- 1
# fd <- "fdisp"
# 
# analyze.fnet(scn.name, run, scn.step, cluster="Dominant", mrc, fd="fdisp", node.importance=F, MM.PC=T)
# 
# 
# 
# 
# scn.name <- "BAU_rcp45"  
# scn.step <- 90
# cluster <- "Dominant"; mrc <- NA; run <- 1; fd <- "fdisp"
# for(cluster in c("Coord", "Dominant")){
#   for(mrc in c(NA, 1:5)){  
#     for(run in 1:3){
#       for(fd in c("fdisp", "fdiv", "fdivw"))
#         analyze.fnet(scn.name, run, scn.step, cluster, mrc, fd, node.importance=F, MM.PC=F)
#     }    
#   }
# }


