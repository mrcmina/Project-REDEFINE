#########################################################################################
## 26
## Plot network properties of the different scenarios
## 
rm(list=ls());gc()
library(tidyverse); library(gridExtra)
library(scales) ; library(ggsci); library(ggpubr)  
cbPalette <- c("#999999", "#56B4E9", "#E69F00",  "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #colorblind-friendly palette
climpal3 <- c("#3CB371", "#8968CD", "#EE9572") # show_col(climpal3)

# show_col(cbPalette)

# Load dataframe with summarized network properties 
prop.df <- read.table("DataOut/SummaryAllNetProp.txt", header=T)
prop.df$mrc[is.na(prop.df$mrc)] <- 0  #replace MRC=NA (whole landscape) with 0

prop.df <- prop.df %>% dplyr::select(.,  -net, -mgmt, -dist, -fd, -Q, -EC) %>% 
             mutate(mgmt_dist=factor(mgmt_dist, levels=c("NONE_NONE","BAU_NONE","BAU_HIGH")), 
             clim=factor(clim, levels=c("CURRENT","RCP45","RCP85")))   

prop.df.means <-  prop.df %>%  group_by(clim, mgmt_dist, mrc, year) %>% 
                      summarise_all(list(mean=mean))

  
# functional diversity 
ggplot(prop.df.means, aes(x=year, y=fdiv_mean,  group=interaction(clim, mgmt_dist), colour=clim, fill=clim)) + 
  geom_line(size=1) +  #geom_ribbon(aes(ymin=fdiv_mean-fdiv_sd, ymax=fdiv_mean+fdiv_sd), alpha = 0.3) +
  facet_grid(mgmt_dist~mrc) + theme_bw()  + 
  scale_color_manual(values=cbPalette)+ scale_fill_manual(values=cbPalette)  +
  labs(y="fdiv", title=paste0("functional diversity ") )

# functional dispersion
ggplot(prop.df.means, aes(x=year, y=fdisp_mean,  group=interaction(clim, mgmt_dist), colour=clim, fill=clim)) + 
  geom_line(size=1) +  #geom_ribbon(aes(ymin=fdisp_mean-fdisp_sd, ymax=fdisp_mean+fdisp_sd), alpha = 0.3) +
  facet_grid(mgmt_dist~mrc) + theme_bw()  + 
  scale_color_manual(values=cbPalette)+ scale_fill_manual(values=cbPalette)  +
  labs(y="fdisp", title=paste0("functional dispersion") )

# functional redundancy 
ggplot(prop.df.means, aes(x=year, y=func.redund_mean,  group=interaction(clim, mgmt_dist), colour=clim, fill=clim)) + 
  geom_line(size=1) +  #geom_ribbon(aes(ymin=func.redund_mean-func.redund_sd, ymax=func.redund_mean+func.redund_sd), alpha = 0.3) +
  facet_grid(mgmt_dist~mrc) + theme_bw()  + 
  scale_color_manual(values=cbPalette)+ scale_fill_manual(values=cbPalette)  +
  labs(y="funct.redund", title=paste0("funcional redundancy")) 

# functional connectivity
ggplot(prop.df.means, aes(x=year, y=PC_mean,  group=interaction(clim, mgmt_dist), colour=clim, fill=clim)) + 
  geom_line(size=1) +  #geom_ribbon(aes(ymin=PC_mean-PC_sd, ymax=PC_mean+PC_sd), alpha = 0.3) +
  facet_grid(mgmt_dist~mrc, scales="free_y") + theme_bw()  + 
  scale_color_manual(values=cbPalette)+ scale_fill_manual(values=cbPalette)  +
  labs(y="funct.conn", title=paste0("functional connectivity"))


## Rescale network properties into percentages and plot them 
mgmt.dist.scens <- c("NONE_NONE","BAU_NONE","BAU_HIGH")
clim.scens <- c("CURRENT","RCP45","RCP85")

mrcs.codes <- seq(0,5,1)
mrcs.name <- c("Centre-du-Quebec", "NicoletYamaska", "Drummond", "Becancour", "LErable","Arthabaska")

runs <- seq(1,5,1)

prop.df.rd <- data.frame()

for (m in 1:length(mgmt.dist.scens)){
  for (c in 1:length(clim.scens)){
    for(y in 1:length(mrcs.codes)){
      for(r in 1:length(runs)){


  mgmt.dist <- mgmt.dist.scens[m]  #"NONE_NONE"
  clim.scen <- clim.scens[c]   #"CURRENT"
  mrc.code <-  mrcs.codes[y] #0
  mrc.name <-  mrcs.name[y]  #"Centre-du-Quebec"
  
  run <- runs[r]

  prop.df.resc <- prop.df %>% filter(mrc==mrc.code, mgmt_dist==mgmt.dist, clim==clim.scen, run==run)   %>%
                    group_by(clim, mgmt_dist, mrc) %>%
                    mutate(fdiv_RD= fdiv/ fdiv[year==2010]-1,
                           fdisp_RD= fdisp/ fdisp[year==2010]-1,
                           func.redund_RD= func.redund/ func.redund[year==2010]-1,
                           PC_RD= PC/ PC[year==2010]-1 ) %>%
                    mutate(fdiv_RD= fdiv_RD*100, 
                           fdisp_RD= fdisp_RD*100, 
                           func.redund_RD= func.redund_RD*100, 
                           PC_RD= PC_RD*100 ) 
    
  
  prop.df.rd <- bind_rows(prop.df.rd, prop.df.resc  )
  
   }
  }
 }
}

prop.df.rd <- prop.df.rd %>% select(.,-fdiv, -fdisp, -func.redund, -PC) %>%
                 group_by(clim, mgmt_dist, year, mrc) %>% 
                 summarise_all(list(mean=mean, sd=sd)) %>% 
                 dplyr::select(., -run_mean,  -run_sd)  


## Rename the scenarios for plotting
prop.df.rd$mgmt_dist <- recode(prop.df.rd$mgmt_dist, NONE_NONE = "CONTROL")
prop.df.rd$mgmt_dist <- recode(prop.df.rd$mgmt_dist, BAU_NONE = "BAU")
prop.df.rd$mgmt_dist <- recode(prop.df.rd$mgmt_dist, BAU_HIGH = "BAU-DIST")

prop.df.rd$clim <- recode(prop.df.rd$clim, CURRENT = "Current")
prop.df.rd$clim <- recode(prop.df.rd$clim, RCP45 = "Moderate")
prop.df.rd$clim <- recode(prop.df.rd$clim, RCP85 = "High")



### PLOT THE THREE INDICATORS AT LANDSCAPE SCALE  
mrc.code <-  0 
mrc.name <-  "Centre-du-Quebec"

prop.df.0 <- prop.df.rd %>% filter(mrc==mrc.code) %>% select(., -mrc) %>%
                rename( climate = clim)

ggplot(prop.df.0, aes(x=year, y=fdiv_RD_mean,  group=interaction(climate, mgmt_dist), fill= climate)) + 
  geom_line() +  geom_ribbon(aes(ymin=fdiv_RD_mean-fdiv_RD_sd, ymax=fdiv_RD_mean+fdiv_RD_sd), alpha = 0.5) +
  facet_grid(~mgmt_dist, scale="free_y") + theme_bw()  + geom_hline(yintercept = 0, colour="gray50", linetype="dashed") +
  scale_color_manual(values=cbPalette)+ scale_fill_manual(values=climpal3)  +
  labs(y="funct.diversity", x="", title= paste0("MRC ", mrc.code," - ",mrc.name)) +  
  theme(strip.text.x = element_text(face="bold"),strip.text.y = element_text( face="bold"),
        strip.background = element_rect(colour="transparent", fill="gray95"))

scaleFUN <- function(x) sprintf("%.1f", x)

p1 <- ggplot(prop.df.0, aes(x=year, y=fdisp_RD_mean,  group=interaction(climate, mgmt_dist), fill= climate, colour=climate)) + 
       geom_line() +  geom_ribbon(aes(ymin=fdisp_RD_mean-fdisp_RD_sd, ymax=fdisp_RD_mean+fdisp_RD_sd), alpha = 0.5) +
       facet_grid(~mgmt_dist, scale="free_y") + theme_pubr()  + geom_hline(yintercept = 0, colour="gray50", linetype="dashed") +
       scale_color_manual(values=climpal3)+ scale_fill_manual(values=climpal3)  +
       labs(y="FDiv", x=element_blank()) +  
       theme(strip.text.x = element_text(face="bold"),strip.text.y = element_text( face="bold"),
             axis.text.x=element_blank(), strip.background = element_rect(colour="transparent", fill="gray95"),
             legend.position = "none", panel.grid.major = element_line(colour="#f0f0f0"), 
             axis.title = element_text(face = "bold",size = rel(1))) + scale_y_continuous(labels=scaleFUN)

p2 <- ggplot(prop.df.0, aes(x=year, y=func.redund_RD_mean,  group=interaction(climate, mgmt_dist), fill=climate, colour=climate)) + 
        geom_line() +  geom_ribbon(aes(ymin=func.redund_RD_mean-func.redund_RD_sd, ymax=func.redund_RD_mean+func.redund_RD_sd), alpha = 0.3) +
        facet_grid(~mgmt_dist, scale="free_y") + theme_pubr()  + geom_hline(yintercept = 0, colour="gray50", linetype="dashed") +
        scale_color_manual(values=climpal3)+ scale_fill_manual(values=climpal3)  +
        labs(y="FRed", x=element_blank()) +  
        theme(strip.text.x = element_blank(), strip.text.y = element_text( face="bold"),
              axis.text.x=element_blank(), strip.background = element_rect(colour="transparent", fill="gray95"),
              legend.position = "none", panel.grid.major = element_line(colour="#f0f0f0"),
              axis.title = element_text(face = "bold",size = rel(1)))

p3 <- ggplot(prop.df.0, aes(x=year, y=PC_RD_mean,  group=interaction(climate, mgmt_dist), fill=climate, colour=climate)) + 
        geom_line() +  geom_ribbon(aes(ymin=PC_RD_mean-PC_RD_sd, ymax=PC_RD_mean+PC_RD_sd), alpha = 0.5) +
        facet_grid(~mgmt_dist, scale="free_y") + theme_pubr()  + geom_hline(yintercept = 0, colour="gray50", linetype="dashed") +
        scale_color_manual(values=climpal3)+ scale_fill_manual(values=climpal3)  +
        labs(y="FConn", x="year")   +  
        theme(strip.text.x = element_blank(), strip.text.y = element_text( face="bold"),
         strip.background = element_rect(colour="transparent", fill="gray95"), panel.grid.major = element_line(colour="#f0f0f0"),
         legend.position = "bottom", legend.title = element_text(face="italic"), 
         axis.title = element_text(face = "bold",size = rel(1)))
 
        
# dir.out <- "C:/Users/marco/OneDrive - UQAM/Documents/REDEFINE/Manuscripts/2 Mina et al TBD/Figures/"
# jpeg(paste0(dir.out,"test.jpeg"), units="in", width=10, height=7, res=600)
grid.arrange(p1, p2, p3,  nrow = 3) 
# dev.off()

g <- arrangeGrob(p1, p2, p3,  nrow = 3)
ggsave(file=paste0(dir.out,"Fnet.prop.CDQ.pdf"), width=10, height=7, g)

# }


## Plot differences between MRCS
prop.net <- prop.df.rd %>% filter(mrc != 0)   #exclude whole CDQ (for now)
                           
# prop.net$mrc <- recode_factor(prop.net$mrc,`0`="Centre-du-Quebec",
#                               `1`="NicoletYamaska", `2`="Drummond", `3`="Becancour",`4`="LErable", `5`="Arthabaska")

prop.net$mrc <- recode_factor(prop.net$mrc,`0`="Centre-du-Quebec",
                              `1`="N.Yam.", `2`="Drum.", `3`="Becan.",`4`="Erabl.", `5`="Arthab.")


#functional dispersion (line charts)
# ggplot(prop.net, aes(year, fdisp_RD_mean, colour = mrc, shape=mrc)) + geom_line(size=1) +
#   geom_point(size = 1) + labs(y="FDisp", x="year")   +
#   facet_grid(mgmt_dist~clim) + theme_pubr() +  scale_color_nejm() +
#   theme( strip.text.x = element_text( face="bold"),strip.text.y = element_text( face="bold"),
#         strip.background = element_rect(colour="transparent", fill="gray95"), panel.grid.major = element_line(colour="#f0f0f0"),
#         legend.position = "bottom", legend.title = element_text(face="italic"),
#         axis.title = element_text(face = "bold",size = rel(1)))

## RADAR PLOTS
library(ggradar)
climpal3rev <- c("#3CB371", "#8968CD", "#EE9572") # show_col(climpal3rev)

#keep only values at the end of the simulation
prop.net.90 <- prop.net  %>% ungroup() %>% filter(year == 2100) 

##-FDisp - CONTROL SCENARIO           
fdisp.90.control.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, fdisp_RD_mean) %>% 
  filter(mgmt_dist=="CONTROL") %>% select(-mgmt_dist) %>% 
  mutate(clim=factor(clim, levels=c("Current","Moderate","High"))) %>% 
  spread(mrc, fdisp_RD_mean)   %>%  rename(group = clim) 

p1 <- ggradar(fdisp.90.control.radar, centre.y=-15, grid.min=-10, grid.mid=0, grid.max=5, values.radar=c("-10%", "0%", "+5%"),
        group.colours=climpal3rev, legend.position="bottom", axis.label.size = 4, group.point.size=4,
        plot.legend=FALSE, grid.label.size=4 ) 


##-FDisp - BAU SCENARIO
fdisp.90.BAU.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, fdisp_RD_mean) %>% 
             filter(mgmt_dist=="BAU") %>% select(-mgmt_dist) %>% 
             spread(mrc, fdisp_RD_mean)   %>% rename(group = clim) 

p2 <- ggradar(fdisp.90.BAU.radar, centre.y=-15, grid.min=-10, grid.mid=0, grid.max=5, values.radar=c("-10%", "0%", "+5%"),
              group.colours=climpal3rev, legend.position="bottom", axis.label.size = 4, group.point.size=4, 
              plot.legend=FALSE, grid.label.size=4 )
 
##-FDisp - BAU-DIST SCENARIO
fdisp.90.BAUDIST.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, fdisp_RD_mean) %>% 
  filter(mgmt_dist=="BAU-DIST") %>% select(-mgmt_dist) %>% spread(mrc, fdisp_RD_mean)   %>% rename(group = clim)

p3 <- ggradar(fdisp.90.BAUDIST.radar, centre.y=-15, grid.min=-10, grid.mid=0, grid.max=5, values.radar=c("-10%", "0%", "+5%"),
              group.colours=climpal3rev, legend.position="bottom", axis.label.size = 4, group.point.size=4,
              plot.legend=FALSE, grid.label.size=4 )



#functional redundancy 
# ggplot(prop.net, aes(year, func.redund_RD_mean, colour = mrc, shape=mrc)) + geom_line(size=1) +
#   geom_point(size = 1) + labs(y="FRed", x="year")   +
#   facet_grid(mgmt_dist~clim) + theme_pubr() +  scale_color_nejm() +
#   theme( strip.text.x = element_text( face="bold"),strip.text.y = element_text( face="bold"),
#          strip.background = element_rect(colour="transparent", fill="gray95"), panel.grid.major = element_line(colour="#f0f0f0"),
#          legend.position = "bottom", legend.title = element_text(face="italic"),
#          axis.title = element_text(face = "bold",size = rel(1)))

##-FRed - CONTROL SCENARIO           
fred.90.control.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, func.redund_RD_mean) %>% 
  filter(mgmt_dist=="CONTROL") %>% select(-mgmt_dist) %>% 
  spread(mrc, func.redund_RD_mean)   %>%  rename(group = clim) 

p4 <- ggradar(fred.90.control.radar, centre.y=-0.3, grid.min=0, grid.mid=0.4, grid.max=0.8, values.radar=c("0%","+0.4%","+0.8%"),
              group.colours=climpal3, plot.legend=FALSE, legend.position="bottom", group.point.size=4, axis.label.size = 4,
              grid.label.size=4)

##-FRed - BAU SCENARIO
fred.90.BAU.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, func.redund_RD_mean) %>% 
  filter(mgmt_dist=="BAU") %>% select(-mgmt_dist) %>% spread(mrc, func.redund_RD_mean)   %>% rename(group = clim) 

p5 <- ggradar(fred.90.BAU.radar, centre.y=-0.3, grid.min=0, grid.mid=0.4, grid.max=0.8, values.radar=c("0%","+0.4%","+0.8%"),
              group.colours=climpal3, plot.legend=FALSE, legend.position="bottom", group.point.size=4, axis.label.size = 4,
              grid.label.size=4)

##-FRed - BAU-DIST SCENARIO
fred.90.BAUDIST.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, func.redund_RD_mean) %>% 
  filter(mgmt_dist=="BAU-DIST") %>% select(-mgmt_dist) %>% spread(mrc, func.redund_RD_mean)   %>% rename(group = clim)

p6 <- ggradar(fred.90.BAUDIST.radar, centre.y=-0.3, grid.min=0, grid.mid=0.4, grid.max=0.8, values.radar=c("0%","+0.4%","+0.8%"),
              group.colours=climpal3, plot.legend=FALSE, legend.position="bottom", group.point.size=4, axis.label.size = 4, 
              grid.label.size=4)

# connectivity 
# ggplot(prop.net, aes(year, PC_RD_mean, colour = mrc, shape=mrc)) + geom_line(size=1) +
#   geom_point(size = 1.5) + labs(y="FConn", x="year")   +
#   facet_grid(mgmt_dist~clim) + theme_pubr() +  scale_color_nejm() +
#   theme( strip.text.x = element_text( face="bold"),strip.text.y = element_text( face="bold"),
#          strip.background = element_rect(colour="transparent", fill="gray95"), panel.grid.major = element_line(colour="#f0f0f0"),
#          legend.position = "bottom", legend.title = element_text(face="italic"),
#          axis.title = element_text(face = "bold",size = rel(1)))


##-FConn - CONTROL SCENARIO           
fconn.90.control.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, PC_RD_mean) %>% 
  filter(mgmt_dist=="CONTROL") %>% select(-mgmt_dist) %>% 
  spread(mrc, PC_RD_mean)   %>%  rename(group = clim) 

p7 <- ggradar(fconn.90.control.radar, centre.y=-50, grid.min=-20, grid.mid=0, grid.max=20, values.radar=c("-20%","0%","+20%"),
        group.colours=climpal3, plot.legend=FALSE, legend.position="bottom", group.point.size=4, axis.label.size = 4, 
        grid.label.size=4)

##-FConn - BAU SCENARIO
fconn.90.BAU.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, PC_RD_mean) %>% 
  filter(mgmt_dist=="BAU") %>% select(-mgmt_dist) %>% spread(mrc, PC_RD_mean)   %>% rename(group = clim) 

p8 <- ggradar(fconn.90.BAU.radar, centre.y=-50, grid.min=-20, grid.mid=0, grid.max=20, values.radar=c("-20%","0%","+20%"),
        group.colours=climpal3, plot.legend=FALSE, group.point.size=4, axis.label.size = 4, grid.label.size=4)

##-FConn - BAU-DIST SCENARIO
fconn.90.BAUDIST.radar <- prop.net.90 %>% select(clim, mgmt_dist, mrc, PC_RD_mean) %>% 
  filter(mgmt_dist=="BAU-DIST") %>% select(-mgmt_dist) %>% spread(mrc, PC_RD_mean)   %>% rename(group = clim)

p9 <- ggradar(fconn.90.BAUDIST.radar, centre.y=-50, grid.min=-20, grid.mid=0, grid.max=20, values.radar=c("-20%","0%","+20%"),
        group.colours=climpal3, plot.legend=FALSE, group.point.size=4, axis.label.size = 4, grid.label.size=4) 

#save in high resolution for the paper
jpeg(paste0(dir.out,"test.jpg"), units="in", width=10, height=8, res=600)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol=3  ) 
dev.off()

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol=3 )
ggsave(file=paste0(dir.out,"Fnet.prop.RCM.pdf"), width=10, height=8, g)

#legend

jpeg(paste0(dir.out,"test.jpg"), units="in", width=10, height=8, res=600)
ggradar(fconn.90.BAUDIST.radar, centre.y=-90, grid.min=-60, grid.mid=-30, grid.max=0, values.radar=c("-60%","-30%","0%"),
        group.colours=climpal3, group.point.size=4, axis.label.size = 4,grid.label.size=4,
        legend.title="climate", plot.legend=TRUE, legend.position="bottom") 
dev.off()

###functions   
theme_Publication <- function(base_size=14) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size)
      + theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(), 
              axis.line = element_line(colour="black"),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.key.size= unit(0.2, "cm"),
              legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic"),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")
      ))
    
  }
  
scale_fill_Publication <- function(...){
    library(scales)
    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
    
  }
  
scale_colour_Publication <- function(...){
    library(scales)
    discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
    
  }