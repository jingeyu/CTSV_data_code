rm(list = ls())
setwd("~/Desktop/PDAC3")
library(ggplot2)
library(qvalue)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(stringr)

cols <- c(rgb(0, 0, 0, maxColorValue = 255),
          "#FA8072",
          # rgb(241,222,184,maxColorValue = 255),
          rgb(0, 0, 255, maxColorValue = 255),
          rgb(131, 31, 138, maxColorValue = 180),
          rgb(230, 145, 56,maxColorValue = 255),
          
          rgb(50, 183, 195, maxColorValue = 255),
          rgb(105,182,42,maxColorValue = 255),
          rgb(147,169,38,maxColorValue = 255),
          # rgb(233,168,136,maxColorValue = 255),
          rgb(138,146,252, maxColorValue = 255),
          rgb(105, 200, 63,maxColorValue = 255))

thre.alpha <- 0.05

###############################################################
######################### Figure 1 ############################
###############################################################
load(paste0("/Users/jinge.yu/Desktop/PDAC/PDAC-data/Real_A_author_st1.RData"))
mean_y <- apply(Y, 1,mean)
variance_y <- apply(Y, 1, var)
#--------------(a) mean vs variance --------------
png(filename = "figures/mean_var.png", width = 700, height = 500)
par(mai = c(1.2,1.4,0.8,0.5))
plot(mean_y, variance_y,
     xlab = "Mean", ylab = "Variance",
     cex.axis = 1.5,
     cex.lab = 2,
     mgp = c(4,1,0),
     bty = "l", # remove the upper and right border
     col = "#2f72b7")
abline(a = 0, b = 1, cex = 2,col = "red")
dev.off()

#-------- (b) histogram of zero inflation -------
zi_y <- colMeans(Y == 0)
png(filename = "figures/zero_hist.png", width = 700, height = 500)
par(mai = c(1.2,1.4,0.8,0.5))
h <- hist(zi_y, xlim = c(0,1),freq = T,
          breaks = seq(0,1,0.05),
)
h$density <- h$counts/sum(h$counts)
plot(h,freq = FALSE,col = "lightgray",
     axes=FALSE,
     cex.lab = 2,
     xlab = "Zero-inflation rate",
     ylab = "Proportion",
     main = "",
     mgp = c(4,1,0),
)
axis(at = seq(0,1,0.1),
     labels = seq(0,1,0.1),
     side = 1,cex.axis = 1.5)
axis(at = seq(0,0.3,0.1),side = 2 ,
     labels = c(0,paste0(c(10,20,30),"%")),
     cex.axis = 1.5)
dev.off()


###############################################################
######################### Figure 5 ############################
###############################################################
#------ Venn plot ------
pdac <- "A"
st <- 1
load(paste0("/Users/jinge.yu/Desktop/PDAC/PDAC-data/Real_",pdac,"_author_st",st,".RData"))
load(file = paste0("real_result_mac/SE_gene_",pdac,"-",st,"_",100*thre.alpha,".RData"))

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
V <- venn.diagram(x = list( "CTSV
(61)" = SE_zinb2,
"SPARK
(958)" = SE_spark,
"BOOST-GP
(814)" = SE_boost,
"SPARK-X
(1800)" = SE_sparkX,
"SpatialDE
(937)" = SE_spatial

),
height = 2200, 
width = 2200,
resolution = 300,
alpha = rep(0.5,5),
col = cols[c(1,3,4,2,5)],
# cex = c(1,2,3,4,5),
# fontface = 5,
cat.dist=c(0.2, 0.24, 0.23,0.23,0.23),
# cat.pos=c(200, 60, 180,120,60),
cat.cex=c(1.0,1.0,1.0,1.0,1.0),
margin = 0.07,
# filename = NULL
filename = paste0("figures/Venn_A1.png")
)


###############################################################
######################### Figure 6 ############################
###############################################################
#-------- gene spatial expression patetern ------------
## anscombe variance stabilizing transformation: NB
var_stabilize <- function(x, sv = 1) {
    varx = apply(x, 1, var)
    meanx = apply(x, 1, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = sv)))
    return(log(x + 1/(2 * phi)))
}# end func

relative_func <- function(expres) {
    maxd = max(expres) - min(expres)
    rexpr = (expres - min(expres))/maxd
    return(rexpr)
}# end func

plot_pattern <- function(y,point_size,spot.coor,filename){
    ex <- relative_func(y)
    pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
    expl.info <- data.frame(spot.coor, ex)
    colnames(expl.info) <- c("x","y","Y")
    rownames(expl.info) <- paste0("spot_",1:ncol(Y))
    gg <- ggplot(expl.info, aes(x = x, y = y, col = Y))
    pl <- gg + geom_point(size = point_size) +
        scale_color_distiller(palette = "Spectral")+
        theme_bw()+
        theme(legend.position = "none",
              legend.title = element_text(size = 0),
              legend.text=element_text(size=15),
              axis.title.x=element_text(size=40,face = "bold"),
              axis.title.y=element_text(size=40,face = "bold"),
              axis.text.x = element_text(size = 35,face = "bold"),
              axis.text.y = element_text(size = 35,face = "bold")
        ) + labs(x = expression(S[1]), y = expression(S[2])) 
    ggsave(filename, pl, width = 11, height = 10)
}

genes <- c("ATXN2L","MED16","AC073896.4","CLPS","CRP","COL6A2")
var_y <- var_stabilize(Y)
for(g in genes){
    plot_pattern(var_y[g,],8, spot.coor, filename = paste0("figures/A1_",g,".png"))
}