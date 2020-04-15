#########################################################################################
## 25 
## Write all properties of a set of scenarios in a single table
## Then plot time series plots with the network properties of the different scenarios
## 
rm(list=ls());gc()
library(dplyr)
library(ggplot2)
library(gridExtra)
#setwd("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork")
 
cluster <- "Dominant"
fd <- "fdisp"

scenarios <- c(paste0(c("BAU", "BAUDistHigh", "NM"), "_current"),
               paste0(c("BAU", "BAUDistHigh", "NM"), "_CanESM2_rcp45"),
               paste0(c("BAU", "BAUDistHigh", "NM"), "_CanESM2_rcp85"))


for(scn.name in scenarios){
  
  ## Identifying the scenario characteristics
  clim <- ifelse(length(grep("rcp45", scn.name)>0), "RCP45", ifelse(length(grep("rcp85", scn.name)>0), "RCP85", "CURRENT"))
  mgmt <- ifelse(length(grep("BAU", scn.name)>0), "BAU", "NONE")
  dist <- ifelse(length(grep("DistHigh", scn.name)>0), "HIGH", "NONE")
  
  ## For each MRC and run, read the network properties of each time step
  for(mrc in c(NA, 1:5)){
    for(run in 1:5){
      for(y in seq(0, 90, 10)){
        
        if(is.na(mrc))
          load(file=paste0("DataScn/NetProp_", scn.name, "_run", run, "_step", y, "_", cluster, "_", fd, ".Rdata"))
        else
          load(file=paste0("DataScn/NetProp_", scn.name, "_run", run, "_step", y, "_", cluster, "_mrc", mrc, "_", fd, ".Rdata"))
        
        
        if(exists("prop"))
          prop <- rbind(prop,
                        data.frame(net=scn.name, clim=clim, mgmt=mgmt, dist=dist, mgmt_dist=paste0(mgmt,"_",dist)  ,
                                   run=run, mrc=mrc, fd=fd, year=y+2010,fdiv=net.prop$net.fdiv, fdisp=net.prop$net.fdisp, 
                                   func.redund=net.prop$net.fredund, 
                                   Q= if(is.null(net.prop$met.modularity))
                                     net.prop$met.modularity=NA
                                   else 
                                     net.prop$met.modularity=net.prop$met.modularity ,
                                   PC=net.prop$PC, EC=net.prop$EC) )
        else

        prop <- data.frame(net=scn.name, clim=clim, mgmt=mgmt, dist=dist, 
                           mgmt_dist=paste0(mgmt,"_",dist)  ,
                           run=run, mrc=mrc, fd=fd, year=y+2010,
                           fdiv=net.prop$net.fdiv, fdisp=net.prop$net.fdisp, func.redund=net.prop$net.fredund, 
                           Q= if(is.null(net.prop$met.modularity))
                             net.prop$met.modularity=NA
                              else 
                                net.prop$met.modularity=net.prop$met.modularity ,
                           PC=net.prop$PC, EC=net.prop$EC)
        
        print(paste0(scn.name))
        
      }
    }
  }
}
write.table(prop, "DataOut/SummaryAllNetProp.txt", quote=F, row.names = F)



## A. Plot ALL SCENARIOS per MRC
for(mrc in c(0:5)){
  tiff(paste0("DataOut/PlotNetProp_mrc", mrc, ".tiff"), width=1200, height=600)
  p1 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=fdiv)) + 
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +
        geom_point(color="black") + ylab("functional diversity") 
  p2 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=fdisp)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional dispersion") 
  p3 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=func.redund)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional redundancy") 
  p4 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=PC)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional connectivity")
  p5 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=Q)) + 
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("modularity")
  grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
  dev.off()
}


## B. Plot ALL SCNEARIOS per MRC, Y-AXIS FULL RANGE
for(mrc in c(0:5)){
  tiff(paste0("DataOut/PlotNetPropFullRange_mrc", mrc, ".tiff"), width=1200, height=600)
  p1 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=fdiv)) + 
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional diversity") + ylim(c(1,6))
  p2 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=fdisp)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional dispersion")  + ylim(c(0,0.4))
  p3 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=func.redund)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional redundancy") + ylim(c(0,1)) 
  p4 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=PC)) +  
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("functional connectivity") + ylim(c(0,1))
  p5 <- ggplot(data=filter(prop, mrc==mrc, cluster==cluster, fd==fd), aes(x=year, y=Q)) + 
        geom_line(aes(colour=mgmt_dist, linetype=clim), size=1) +  
        geom_point(color="black") + ylab("modularity") + ylim(c(0,1))
  grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
  dev.off()
}


## Plot indexes for ONE scenario per MRC
scn.name <- "BAU_CanESM2_rcp85"
  
tiff(paste0("DataOut/PlotNetProp_", scn.name, ".tiff"), width=1000, height=600)
p1 <- ggplot(data=filter(prop,fd=="fdisp"), aes(x=year, y=fdiv, group=net)) +  geom_line(color=mrc) +  
      geom_point(color="black") + ylab("functional diversity") + ggtitle(scn.name) + theme(plot.title = element_text(size=18))
p2 <- ggplot(data=filter(prop,fd=="fdisp"), aes(x=year, y=fdisp, group=net)) +  geom_line(color=mrc) + 
      geom_point(color="black") + ylab("functional dispersion") + ggtitle("") + theme(plot.title = element_text(size=18))
p3 <- ggplot(data=filter(prop,fd=="fdisp"), aes(x=year, y=func.redund, group=net)) +  geom_line(color=mrc) + 
      geom_point(color="black")  + ylab("functional redundancy") + ggtitle("") + theme(plot.title = element_text(size=18))
p4 <- ggplot(data=filter(prop,fd=="fdisp"), aes(x=year, y=PC, group=net)) +  geom_line(color=mrc) +  
      geom_point(color="black")  + ylab("functional connectivity")
p5 <- ggplot(data=filter(prop,fd=="fdisp"), aes(x=year, y=Q, group=net)) +  geom_line(color=mrc) +  
      geom_point(color="black")  + ylab("modularity")
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
dev.off()




