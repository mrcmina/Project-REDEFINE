#########################################################################################
## 23 
## Four options to plot a functional network
## I.   node size as patch area-----node color as fdiv------link width as pctg.fdiv.w
## II.  node size as patch area-----node color as fdisp-----link width as pctg.fdisp
## III. node size as patch area-----node color as fdiv------link color as pctg.fdiv.w
## IV.  node size as patch area-----node color as fdisp-----link color as pctg.fdiv.w
#########################################################################################

rm(list=ls())
library(igraph); library(viridis)
#setwd("C:/WORK_Offline/REDEFINE_offline/FunctionalNetwork")

## Parameters of the functional network to be plotted

mgmt.scenarios <-  c("NM", "BAU",  "BAUDistHigh")

scenarios <- c(paste0(mgmt.scenarios,"_", "current"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp45"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp85"))


  
  for(scn.name in scenarios){
    for(mrc in c(NA, 1:5)){ 
        for(scn.step in seq(0,90,10)) {


step <- scn.step
run <- 1
cluster <- "Dominant"


## Load it
if(is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
if(!is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".Rdata"))


## The spatial layout: To have the nodes positioned at their centroid, 
## create a normalized matrix with the coordinates
net.coord <- matrix(c(V(func.net)$xm, V(func.net)$ym), nc=2)
layout <- net.coord/max(net.coord)
layout.plot <- net.coord/max(net.coord)


## 4 Colours for fd of nodes
node.colors <- viridis(4)  #show_col(node.colors)
link.colors <- grey.colors(4)[4:1]  #show_col(link.colors)


## I. Plot the network
## Vertex's size as patch's area (in the landscape)
## Vertex's color as functional diversity (4 classes)
## Edge's width as link's weighted pctg func.diversity
## Edge's color constant
# node.size <- (V(func.net)$size)^0.25
# node.score <- ifelse(V(func.net)$fdiv<=2.25,1,ifelse(V(func.net)$fdiv<=3.5,2,ifelse(V(func.net)$fdiv<=4.75,3,4)))
# link.width <- E(func.net)$pctg.fdiv.w*4
# 
# tiff(paste0("DataScn/PlotFnet/PlotNet_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, "_size_fdiv_pctgfdivw.tiff"), width = 800, height = 800)
# plot(func.net, layout=layout, margin=c(0.1,0.1,0.1,0.1),
#      vertex.label=NA, vertex.color=node.colors[node.score],
#      vertex.frame.color=NA, vertex.size=node.size,
#      edge.arrow.size=0, edge.width=link.width, edge.color="grey80")
# legend(-1.2, 1, title="Functional Diversity",
#        legend=c("1 - 2.25", "2.25 - 3.50", "3.50 - 4.75", "4.75 - 6"),
#        pt.cex=2.5, col='black', pch=21, pt.bg=node.colors, bty="n")
# legend(-1.2, -0.3, title="Patch area",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
# legend(-1.2, -0.7, title="Link weight",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=15, pt.bg = "grey", bty="n")
# dev.off()


## II. Plot the network: 
## Vertex's size as patch's area (in the landscape)
## Vertex's color as functional dispersion (4 classes)
## Edge's width as pctg.fdisp
## Edge's color constant
node.size <- (V(func.net)$size)^0.25
node.score <- ifelse(V(func.net)$fdisp<=0.1,1,ifelse(V(func.net)$fdisp<=0.2,2,ifelse(V(func.net)$fdisp<=0.3,3,4)))
link.size <- E(func.net)$pctg.fdisp*4


# tiff(paste0("DataScn/PlotFnet/PlotNet_", scn.name, "_run", run, "_step", step, "_", cluster,  "_mrc", mrc, "_size_fdisp_pctgfidsp.tiff"), width = 800, height = 800)

plot(func.net, layout=layout, margin=c(0.1,0.1,0.1,0.1),
     vertex.label=NA, vertex.color=node.colors[node.score], 
     vertex.frame.color=NA, vertex.size=node.size, 
     edge.arrow.size=0,  edge.width=link.size, edge.color="grey80",
     main=paste0(scn.name," - year ",step+2010)  )

legend(-1.2, 1, title="Functional Dispersion",
       legend=c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.37 - 0.4"),
       pt.cex=2.5, col='black', pch=21, pt.bg=node.colors, bty="n")          
legend(-1.2, -0.3, title="Patch area",
       legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
legend(-1.2, -0.7, title="Link weight",
       legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=15, pt.bg = "grey", bty="n")

# dev.off()


      

        } #loop step
    } #loop mrc
  } #loop scenarios

## Plot using ggnetwork
# library(ggnetwork); library(intergraph)
# func.net.df <- ggnetwork(func.net)
# func.net.df$node.score <- c(node.score, rep(NA, nrow(func.net.df)-length(node.score) ) )
# func.net.df$link.size <- func.net.df$pctg.fdisp*4
# 
# 
# ggplot(func.net.df, aes(x = xm, y = ym, xend = xend, yend = yend), layout=layout.plot) +
#   geom_edges(aes(size=link.size),
#               color = "grey70") +
#   geom_nodes(aes(color = fdisp, size = size)) +
#   theme_blank()



## III. Plot the network: 
## Vertex's size as patch's area (in the landscape)
## Vertex's color as functional diversity (4 classes)
## Edge's width constant
## Edge's color as pctg.fdiv.w
# node.size <- (V(func.net)$size)^0.25
# node.score <- ifelse(V(func.net)$fdiv<=2.25,1,ifelse(V(func.net)$fdiv<=3.5,2,ifelse(V(func.net)$fdiv<=4.75,3,4)))
# q <- quantile(E(func.net)$pctg.fdiv.w)
# link.score <- ifelse(E(func.net)$pctg.fdiv.w<=q[2], 1, ifelse(E(func.net)$pctg.fdiv.w<=q[3], 2, ifelse(E(func.net)$pctg.fdiv.w<=q[3], 3, 4)))
# tiff(paste0("DataScn/PlotFnet/PlotNet_", scn.name, "_run", run, "_step", step, "_", cluster,  "_mrc", mrc, "_size_fdiv_col.pctgfdivw.tiff"), width = 800, height = 800)
# plot(func.net, layout=layout, margin=c(0.1,0.1,0.1,0.1),
#      vertex.label=NA, vertex.color=node.colors[node.score],
#      vertex.frame.color=NA, vertex.size=node.size,
#      edge.arrow.size=0,  edge.color=link.colors[link.score])
# legend(-1.2, 1, title="Functional Diversity",
#        legend=c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4"),
#        pt.cex=2.5, col='black', pch=21, pt.bg=node.colors, bty="n")
# legend(-1.2, -0.3, title="Patch area",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
# legend(-1.2, -0.7, title="Link weight",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col=link.colors, pch=15, pt.bg = "grey", bty="n")
# dev.off()


## IV. Plot the network: 
## Vertex's size as patch's area (in the landscape)
## Vertex's color as functional dispersion (4 classes)
## Edge's width constant
## Edge's color as pctg.fdisp
# node.size <- (V(func.net)$size)^0.25
# node.score <- ifelse(V(func.net)$fdisp<=0.1,1,ifelse(V(func.net)$fdisp<=0.2,2,ifelse(V(func.net)$fdisp<=0.3,3,4)))
# q <- quantile(E(func.net)$pctg.fdisp)
# link.score <- ifelse(E(func.net)$pctg.fdisp<=q[2], 1, ifelse(E(func.net)$pctg.fdisp<=q[3], 2, ifelse(E(func.net)$pctg.fdisp<=q[3], 3, 4)))
# tiff(paste0("DataScn/PlotFnet/PlotNet_", scn.name, "_run", run, "_step", step, "_", cluster,  "_mrc", mrc, "_size_fdisp_col.pctgfdisp.tiff"), width = 800, height = 800)
# plot(func.net, layout=layout, margin=c(0.1,0.1,0.1,0.1),
#      vertex.label=NA, vertex.color=node.colors[node.score], 
#      vertex.frame.color=NA, vertex.size=node.size,
#      edge.arrow.size=0,  edge.color=link.colors[link.score])
# legend(-1.2, 1, title="Functional Dispersion",
#        legend=c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4"),
#        pt.cex=2.5, col='black', pch=21, pt.bg=node.colors, bty="n")           
# legend(-1.2, -0.3, title="Patch area",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
# legend(-1.2, -0.7, title="Link weight",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col=link.colors, pch=15, pt.bg = "grey", bty="n")
# dev.off()


### PLOT FNET for paper 

## Parameters of the functional network to be plotted

mgmt.scenarios <-  c("NM", "BAU",  "BAUDistLow","BAUDistHigh")

scenarios <- c(paste0(mgmt.scenarios,"_", "current"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp45"),
               paste0(mgmt.scenarios,"_", "CanESM2_rcp85"))
    
      
scn.name <- "BAUDistHigh_CanESM2_rcp85"
step <- 0
run <- 1
cluster <- "Dominant"
mrc <- NA
      
      
## Load it
if(is.na(mrc))
      load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
if(!is.na(mrc))
      load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".Rdata"))
    

## The spatial layout: To have the nodes positioned at their centroid, 
## create a normalized matrix with the coordinates
net.coord <- matrix(c(V(func.net)$xm, V(func.net)$ym), nc=2)
layout <- net.coord/max(net.coord)
layout.plot <- net.coord/max(net.coord)
      
      
## 4 Colours for fd of nodes
node.colors <- viridis(4)  #show_col(node.colors)
link.colors <- grey.colors(4)[4:1]  #show_col(link.colors)

## II. Plot the network: 
## Vertex's size as patch's area (in the landscape); Vertex's color as functional dispersion (4 classes)
## Edge's width as pctg.fdisp ; Edge's color constant
node.size <- (V(func.net)$size)^0.25
node.score <- ifelse(V(func.net)$fdisp<=0.1,1,ifelse(V(func.net)$fdisp<=0.2,2,ifelse(V(func.net)$fdisp<=0.3,3,4)))
link.size <- E(func.net)$pctg.fdisp*4
      
      
# tiff(paste0("DataScn/PlotFnet/PlotNet_", scn.name, "_run", run, "_step", step, "_", cluster,  "_mrc", mrc, "_size_fdisp_pctgfidsp.tiff"), width = 800, height = 800)
 
plot(func.net, layout=layout, margin=c(rep(0, 4) ),
           vertex.label=NA, vertex.color=node.colors[node.score],
           vertex.frame.color=NA, vertex.size=node.size,
           edge.arrow.size=0,  edge.width=link.size, edge.color="gray80")

legend(-1.2, 1, title="Functional Dispersion",
             legend=c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.37 - 0.4"),
             pt.cex=2.5, col='black', pch=21, pt.bg=node.colors, bty="n")          
legend(-1.2, -0.3, title="Patch area",
             legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
legend(-1.2, -0.7, title="Link weight",
             legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=15, pt.bg = "grey", bty="n")
      
# dev.off()



## Plot using ggraph or ggnetwork
# library(ggraph); library(ggnetwork)
# 
# cols_f <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))
# 
# ggraph(func.net, layout=layout.plot) + 
#   geom_node_point(aes(colour = as.factor(node.colors[node.score]), size=node.size  )) +
#   scale_color_manual(
#     limits = as.factor(node.colors[node.score]),
#     values = cols_f(nrow(layout.plot))
#   ) +
#   coord_equal()+ 
#   geom_edge_link(colour = "grey80") + 
#   theme_graph()




### Load two networks and compute difference between them 
library(tidygraph)

## Parameters of the functional network to be plotted
run <- 1 ; cluster <- "Dominant" ; mrc <- 4  #common betweeen the two networks 


 
# for(scn.step2 in c(0,90)) {
  
scn.step2 <- 90

#First network 
# scn.name <- "NM_current"
# scn.title <- "CONTROL - Current"
  
  # scn.name <- "NM_CanESM2_rcp45"
  # scn.title <- "CONTROL - Moderate"
 
  # scn.name <- "NM_CanESM2_rcp85"
  # scn.title <- "CONTROL - High"

# scn.name <- "BAU_CanESM2_rcp45"
# scn.title <- "BAU - MODERATE"
  
  # scn.name <- "BAU_current"
  # scn.title <- "BAU - CURRENT"
  
  # scn.name <- "BAU_CanESM2_rcp85"
  # scn.title <- "BAU - HIGH"

scn.name <- "BAUDistHigh_CanESM2_rcp85"
scn.title <- "BAU-DIST - HIGH"
  
  # scn.name <- "BAUDistHigh_current"
  # scn.title <- "BAU-DIST - CURRENT"
  
  # scn.name <- "BAUDistHigh_CanESM2_rcp45"
  # scn.title <- "BAU-DIST - MODERATE"
  

step <- 0
step2 <- scn.step2

if(is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
if(!is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".Rdata"))

func.net1 <- func.net #save it as object 
# func.net1.df <- as_long_data_frame(func.net1)

func.net1_tidy <- as_tbl_graph(func.net1)
net1 <- as.data.frame(func.net1_tidy)
fdisp.net1 <- net1$fdisp                 #save the nodes' fdisp of the first network as object 
# hist(fdisp.net1)

## Plot the first network 
net.coord <- matrix(c(V(func.net1)$xm, V(func.net1)$ym), nc=2)
layout <- net.coord/max(net.coord)
layout.plot <- net.coord/max(net.coord)
node.colors <- viridis(4)  #show_col(node.colors)
link.colors <- grey.colors(4)[4:1]  #show_col(link.colors)
node.size <- (V(func.net1)$size)^0.25
node.score <- ifelse(V(func.net1)$fdisp<=0.1,1,ifelse(V(func.net1)$fdisp<=0.2,2,ifelse(V(func.net1)$fdisp<=0.3,3,4)))
link.size <- E(func.net1)$pctg.fdisp*4
link.width <- E(func.net)$pctg.fdiv.w*4

jpeg(paste0("DataScn/PlotFnet/netdiff/supplementary_scen/PlotNet_",scn.name, "_step", step2, "_mrc", mrc, ".jpg"), width = 30, height = 30, unit="cm", res=300 )
plot(func.net1, layout=layout, margin=c(rep(0, 4) ), vertex.label=NA, vertex.color=node.colors[node.score], 
             vertex.frame.color=NA, vertex.size=node.size, edge.arrow.size=0.4,  edge.color="gray85",
             edge.width=link.width) #"#e2a4a4"


## SECOND network 
step <- step2

if(is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, ".Rdata"))
if(!is.na(mrc))
  load(file=paste0("DataScn/FuncNetFDISP_", scn.name, "_run", run, "_step", step, "_", cluster, "_mrc", mrc, ".Rdata"))

func.net2 <- func.net  #save it as object 

func.net2_tidy <- as_tbl_graph(func.net2) 
net2 <- as.data.frame(func.net2_tidy)
fdisp.net2 <- net2$fdisp                   #save nodes' fdisp of the second network as object
#hist(fdisp.net2)
fdisp.net12 <- fdisp.net1 - fdisp.net2     #calculate difference in fdisp between the two networks
# hist(fdisp.net12, breaks=20)

#change the nodes FDisp of second network with the difference between the two FDisp
func.net.diff_tidy <- func.net2_tidy %>% activate(nodes) %>% mutate(fdisp = fdisp.net1 - fdisp)  
func.net.diff <- as.igraph(func.net.diff_tidy) #covert back into igraph for plotting

## The spatial layout: To have the nodes positioned at their centroid, create a normalized matrix with the coordinates
net.coord <- matrix(c(V(func.net.diff)$xm, V(func.net.diff)$ym), nc=2)
layout.plot <- net.coord/max(net.coord)

## Colours for FDiv of nodes
node.colors <- c("firebrick3", "gray65", "#35B779FF")  #show_col(node.colors)
# node.colors.y <- c("firebrick3", "#EEE8AA", "#00CD66")  #show_col(node.colors.y)
link.colors <- grey.colors(4)[4:1]  #show_col(link.colors)


## II. Plot the network: 
## Vertex's size as patch's area (in the landscape); Vertex's color as functional dispersion (4 classes)
## Edge's width as pctg.fdisp ; Edge's color constant
node.size <- (V(func.net.diff)$size)^0.25
node.score <- ifelse(V(func.net.diff)$fdisp <  (0) ,  1,
              ifelse(V(func.net.diff)$fdisp <= 0.1,  2,
              ifelse(V(func.net.diff)$fdisp >  0.1,  3, 4)))


link.size <- E(func.net.diff)$pctg.fdisp*4
link.width <- E(func.net)$pctg.fdiv.w*4


par(new=TRUE) #to overplot 
plot(func.net.diff, layout=layout, margin=c(rep(0, 4) ),
        vertex.label=NA, vertex.color=node.colors[node.score], 
        vertex.frame.color=NA, vertex.size=node.size, edge.width=link.width,
        edge.arrow.size=0.4,   edge.color="grey40", main=paste0(scn.title," - year ",step+2010) )

 dev.off()

# }


# 
# legend("bottomright", legend=c("< 0", "0-0.1)", "> 0.1", "links 2010", "links 2100"),
#        pt.cex=2, col=c('black','black','black',"grey85","grey40"), pch=c(21,21,21,-8594, -8594), 
#        pt.bg=c("firebrick3","gray65","#35B779FF","grey85","grey40"), bty="n")  
# 
# 
# 
# plot.new()
# legend("center", 
#        legend=c("< 0", "0-0.1", ">0.1"),
#        pt.cex=2, col=c('black','black','black'), pch=c(21,21,21), 
#        pt.bg=c("firebrick3","gray65","#35B779FF"), bty="n", horiz=TRUE)  
# 
# legend("center", 
#        legend=c( " 2010", " 2100"),
#        pt.cex=2, col=c("grey85","grey40"), pch=c(-8594, -8594), 
#        pt.bg=c("grey85","grey40"), bty="n", horiz=TRUE)  


# legend(-1.2, -0.3, title="Patch area",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=21, pt.bg = "black", bty="n")
# legend(-1.2, -0.7, title="Link weight",
#        legend=c("-", "", "", "+"), pt.cex=c(1,1.5,2,3), col="black", pch=15, pt.bg = "grey", bty="n")



# library(visNetwork)
# data <- toVisNetworkData(func.net.diff)
# data$nodes$size <- data$nodes$size^0.25
# visNetwork(nodes = data$nodes, edges = data$edges) %>%
#   visIgraphLayout() %>%
#   visNodes(size = 10)
