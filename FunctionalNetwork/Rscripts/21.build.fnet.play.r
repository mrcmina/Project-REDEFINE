#########################################################################################
## 21
## Run the function that builds the functional network:
## It saves .Rdata with nodes and links of the network, and the network itself
## It also writes .txt needed to analyze connectivity with conefor.exe software
#########################################################################################

rm(list=ls());gc()
# setwd("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork") # setwd("C:/WORK/REDEFINE")
source("Rscripts/21.build.fnet.r")

## Set of final scenarios with replicates 
scn.path <- "F:/LANDIS_PnET/simulations/Outputs/" # the path to LANDIS simulations (e.g., external drive)

mgmt.scenarios <-  c("NM", "BAU", "BAUDistHigh")  

scenarios <- c(paste0(mgmt.scenarios,"_", "CanESM2_rcp45"),  #paste0(mgmt.scenarios,"_", "current"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp85"))


for(scn.name in scenarios){
  for(run in 1:5){
    for(step in seq(0,90,10)){
      for(mrc in c(NA, 1:5))  {
        
       build.fnet(scn.path, scn.name, run, step, cluster="Dominant", mrc)
      
    } #mrc
  } #step
 } #run
} #scen




## Test for one scenario
scn.name="BAUDistHigh_current"
run=2
step=40
cluster="Dominant"
mrc="5"
build.fnet(scn.path, scn.name, run, step, cluster="Dominant", mrc)



# ## Analyze 3 replicates of the same scenarios, the final snapshop, 
# ## under 2 modes or apporaches of clustering to define patches
# scn.path <- "LANDIS/CdQsims/Output/"  
# scn.name <- "BAU_rcp45"  
# build.fnet(scn.path, scn.name, 1, step=90, cluster="Coord", mrc=NA)
# build.fnet(scn.path, scn.name, 2, step=90, cluster="Coord", MRC=NA)
# build.fnet(scn.path, scn.name, 3, step=90, cluster="Coord", MRC=NA)
# build.fnet(scn.path, scn.name, 1, step=90, cluster="Dominant", MRC=NA)
# build.fnet(scn.path, scn.name, 2, step=90, cluster="Dominant", MRC=NA)
# build.fnet(scn.path, scn.name, 3, step=90, cluster="Dominant", MRC=NA)
# 
# 
# ## 07/10/2019
# ## Analyze the functional networks drawn from the MRC division
# scn.path <- "LANDIS/CdQsims/Output/"  
# scn.name <- "BAU_rcp45" 
# for(mrc in c(NA,1:5)){
#   for(run in 1:3){
#     build.fnet(scn.path, scn.name, run, step=90, cluster="Coord", mrc)
#     build.fnet(scn.path, scn.name, run, step=90, cluster="Dominant", mrc)
#   }
# }

