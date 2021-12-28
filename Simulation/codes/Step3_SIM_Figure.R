rm(list = ls())
setwd("~/Desktop/SIM")
library(ggplot2)
library(pheatmap)
library(pROC)

###############################################################
######################### Figure 2 ############################
###############################################################
#--------------(a) Spot region pattern --------------
load("data/spot_filter_ct4.RData")
pal2 <- c(rgb(200, 50, 54, maxColorValue = 200), 
          rgb(224, 190, 56, maxColorValue = 230),
          rgb(50, 183, 195, maxColorValue = 255), 
          rgb(131, 31, 138, maxColorValue = 180))
pal2 <- setNames(pal2, c("1", "2", "3", "4"))

spot.coor <- as.data.frame(spot.coor)
gg <- ggplot(spot.coor, aes(x = x, y = y))
pl <- gg + geom_point(size = 3, 
                      aes(color = as.factor(spot.region))) +
    scale_color_manual(values=pal2) +
    theme(legend.text=element_text(size=15),        
          legend.position = "top", 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          panel.background = element_blank()) +
    guides(color = guide_legend(title = "Region", 
                                title.theme = element_text(size = 10),
                                override.aes = list(size = 7)))

ggsave(paste0("figures/region.png"), pl, width = 5, height = 4)

#-------- (b) SV gene pattern heatmap --------
load("data/Sim_linear.RData")
G <- nrow(Y)
n <- ncol(Y)
K <- nrow(W)
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(2)
#True cell effects
anno_col <- data.frame(Celltype = factor(c("1", "2", "3","4","5","6")))
colnames(gamma_true) <- c("C1", "C2", "C3","C4","C5","C6")
# Notice that the colnames of matrix and anno_col must be the same!
rownames(anno_col) <- colnames(gamma_true)
anno_cor <- list(Celltype = c("1"= rgb(112, 173, 71, maxColorValue = 255),
                              "2" = rgb(47, 85, 151, maxColorValue = 255),
                              "3" =  rgb(237, 125, 49, maxColorValue = 255),
                              "4" = rgb(224, 190, 56, maxColorValue = 230),
                              "5" = rgb(50, 183, 195, maxColorValue = 255), 
                              "6" = rgb(131, 31, 138, maxColorValue = 180)))

# rownames(gamma_true) <- paste0("Gene",1:G)
labels_row <- c("- 1",rep("",48))
for(lab in 1:(1000/50-1)){
    labels_row = c(labels_row, paste0("- ",lab*50),rep("",49))
}
labels_row <- c(labels_row,paste0("-",1000))


p1 <- pheatmap(gamma_true[1:1000,],
               color = colors,
               # breaks = seq(0,1,length.out = length(colors)+1),
               cluster_cols = F, cluster_rows = F,
               show_rownames = T, show_colnames = F,
               legend_breaks = c(0,1),
               annotation_col = anno_col,
               annotation_colors = anno_cor,
               annotation_names_col = F,
               annotation_legend = F,
               # show_rownames = T,
               labels_row = labels_row,
               fontsize_row = 6,
               width = 3, height = 2.5,
               filename = paste0("figures/gamma_true.png")
)

###############################################################
######################### Figure 3 ############################
###############################################################
#---- (a) ROC curves refer to sim_ana_mac.R -----
pattern <- c("linear","gau","cos")
true_se <- rep(0,G)
true_se[which(rowSums(gamma_true) > 0)] <- 1
cols <- c(rgb(0, 0, 0, maxColorValue = 255),
          rgb(255, 0, 0, maxColorValue = 255),
          rgb(0, 0, 255, maxColorValue = 255),
          rgb(131, 31, 138, maxColorValue = 180),
          rgb(230, 145, 56,maxColorValue = 255),
          rgb(50, 183, 195, maxColorValue = 255),
          rgb(105,182,42,maxColorValue = 255),
          rgb(147,169,38,maxColorValue = 255),
          rgb(138,146,252, maxColorValue = 255),
          rgb(105, 200, 63,maxColorValue = 255))
for(pat in pattern){
   P_val_ZINB <- read.csv(paste0("re/CTSV_sim_",pat,"_pval.csv"))
   rownames(P_val_ZINB) <- P_val_ZINB[,1]
   P_val_ZINB <- P_val_ZINB[,-1]
   # Q_ZINB1 <- read.csv(paste0("re/CTSV_sim_",pat,"_qval.csv"))
   # rownames(Q_ZINB1) <- Q_ZINB1[,1]
   # Q_ZINB1 <- Q_ZINB1[,-1]
   P_val_ZINB1<- matrix(NA,G,K)
    for(k in 1:K){
        P_val_ZINB1[,k] <- apply(P_val_ZINB[,c(k,k+K)],1,min)
    }
    rownames(P_val_ZINB) <- rownames(P_val_ZINB)
    p_val_zinb1 <- apply(P_val_ZINB1,1,min)

    trend_se <- read.csv(paste0("re/trend_sim_",pat,".csv"), header = TRUE)
    rownames(trend_se) <- trend_se[,1]
    trend_se <- trend_se[,-1]
    p_val_trend <- apply(trend_se, 1, min)
    # print(sum(p_val_trend < 0.05))
    spark_se <- read.csv(paste0("re/spark_sim_",pat,".csv"), header = TRUE)
    sparkx_se <- read.csv(paste0("re/sparkX_sim_",pat,".csv"), header = TRUE)
    spatial_se <- read.csv(paste0("re/spatial_sim_",pat,".csv"))
    if(length(which(duplicated(spatial_se) == 1)) > 0)
        spatial_se <- spatial_se[-which(duplicated(spatial_se) == 1),]
    spatial_gene <- spatial_se[,1]
    boost_se <- read.csv(paste0("re/BOOST_sim_",pat,".csv"))
    boost_ppi <- boost_se$PPI
    rownames(boost_se) <- paste0("Gene",1:G)
    somde_se <- read.csv(paste0("re/somde_sim_",pat,".csv"))
    somde_gene <- somde_se[,1]
    # print(sum(as.matrix(somde_se)[,3] < 0.05))
    
    #----roc CURVE-----
    p_val_sparkx <- sparkx_se[,2]
    p_val_spatialDE <- spatial_se$pval
    p_val_spark <- spark_se[,2]
    ind3 <- match(rownames(P_val_ZINB), spatial_se$X) 
    p_val_spatialDE <- as.numeric(spatial_se$pval[ind3])
    ind4 <- match(rownames(P_val_ZINB), somde_se[ ,1]) 
    p_val_somde <- as.numeric(somde_se$pval[ind4])
    p_val_trend_emark <- trend_se[,1]
    p_val_trend_markcorr <- trend_se[,2]
    p_val_trend_markvario <- trend_se[,3]
    p_val_trend_vmark <- trend_se[,4]
    
    roc_spark <- roc(true_se~p_val_spark,algorithm = 1,direction = ">")
    roc_sparkx <- roc(true_se~p_val_sparkx,algorithm = 1,direction = ">")
    roc_spatialDE <- roc(true_se~p_val_spatialDE,algorithm = 1,direction = ">")
    roc_zinb1 <- roc(true_se~p_val_zinb1,algorithm = 1,direction = ">")
    roc_boost <- roc(true_se~boost_ppi, algorithm = 1, direction = ">")
    roc_emark <- roc(true_se~p_val_trend_emark,algorithm = 1,direction = ">")    
    roc_mc <- roc(true_se~p_val_trend_markcorr,algorithm = 1,direction = ">")    
    roc_mv <- roc(true_se~p_val_trend_markvario,algorithm = 1,direction = ">")    
    roc_vmark <- roc(true_se~p_val_trend_vmark,algorithm = 1,direction = ">")
    roc_somde <- roc(true_se~p_val_somde,algorithm = 1,direction = ">")
    
    #-----    
    lwd = 2
    png(filename = paste0("figures/ROC_CG_",pat,".png"),
        width = 1100, height = 1000, res = 200)
    par(mai=c(0.9,0.9,0.2,0.5))
    # par(mai=c(0.9,0.9,0.2,2.0))
    xlim0 <- 0.05
    plot(1-roc_zinb1$specificities,roc_zinb1$sensitivities,
         col = cols[1],lwd = lwd,type = "l",
         xlim = c(0,xlim0),ylim = c(0,1), 
         xlab = "FPR", ylab = "TPR",
         mgp=c(3, 1, 0), lty = 1,
         cex.lab = 1.5)
    
    lines(1-roc_sparkx$specificities,roc_sparkx$sensitivities,
          type = "l",col = cols[2],lwd = lwd,
          lty = 2,xlim = c(0,xlim0))
    
    lines(1-roc_spark$specificities,roc_spark$sensitivities,
          type = "l",col = cols[3],lwd = lwd,
          lty = 3,xlim = c(0,xlim0))
    
    lines(1-roc_boost$specificities,roc_boost$sensitivities,
          type = "l",col = cols[4],lwd = lwd,
          lty = 4,xlim = c(0,xlim0))
    
    lines(1-roc_spatialDE$specificities,roc_spatialDE$sensitivities,
          type = "l",col = cols[5],lwd = lwd,
          lty = 5,xlim = c(0,xlim0))
    lines(1-roc_mv$specificities,roc_mv$sensitivities,
          type = "l",col = cols[6],lwd = lwd,
          lty = 6,xlim = c(0,xlim0))
    lines(1-roc_somde$specificities,roc_somde$sensitivities,
          type = "l",col = cols[7],lwd = lwd,
          lty = 7,xlim = c(0,xlim0))
    
    # legend(xlim0+0.005,0.6, c("ZINB","SPARK-X","SPARK","BOOST-GP","SpatialDE",
    #                           "trendsceek",
    #                           "SOMDE"),
    # col = cols[c(1:7)], cex = 1, lwd = 1,lty = c(1:7),
    # xpd = TRUE)
    dev.off()
}

#---- (b) SV gene expression patterns-----
pal3 <- colorRampPalette(c("#fbecd9", "#4f9a8a"))

for(pat in pattern){
    load(paste0("data/Sim_",pat,".RData"))
    S <- t(spot.coor) - colMeans(spot.coor)
    S <- t(S / apply(S, 1, sd))
    g <- 300
    k <- 1
    mu_gk <- eta[g,k] + beta1[g,k] * h1 + beta2[g,k] * h2
    expl.info <- data.frame(S, mu_gk)
    colnames(expl.info) <- c("x","y","mu")
    rownames(expl.info) <- paste0("spot_",1:length(mu_gk))
    gg <- ggplot(expl.info, aes(x = x, y = y, col = mu))
  
    pl <- gg + geom_point(size = 10) +
        scale_color_gradientn(colours = pal3(5))+
        theme_bw()+
        theme(legend.position = "none",
              # legend.text=element_text(size=20, face = 'bold'),
              axis.title.x=element_text(size=45,face = "bold"),
              axis.title.y=element_text(size=45,face = "bold"),
              axis.text.x = element_text(size = 40,face = "bold"),
              axis.text.y = element_text(size = 40,face = "bold")
        ) +
        labs(x = expression(S[1]), y = expression(S[2])) 
    ggsave(paste0("figures/simulation_",pat,".png"), pl, width = 12, height = 8)
    
    
}

###############################################################
######################### Figure 4 ############################
###############################################################
colors_func <- colorRampPalette(c('white', "black"))
colors <- colors_func(100)
for(pat in pattern){
    P_val <- read.csv(file = paste0("re/CTSV_sim_",pat,"_pval.csv"))
    rownames(P_val) <- P_val[,1]
    P_val <- P_val[,-1]

    P_VAL <- matrix(NA,G,K)
    P_sp <- matrix(NA,G,K)
    for(k in 1:K){
        P_VAL[,k] <- apply(P_val[,c(k,k+K)],1,min)
        spx_ct <- read.csv(paste0("re/sparkX_sim_",pat,"_ct_",k,".csv"))
        P_sp[,k] <- spx_ct$combinedPval
    }
    
    logp_ctsv <- -log10(P_VAL)
    logp_sparkx <- -log10(P_sp)
    

    #True cell effects
    anno_col <- data.frame(Celltype = factor(c("1", "2", "3","4","5","6")))
    # Notice that the colnames of matrix and anno_col must be the same!
    colnames(logp_ctsv) <- c("C1", "C2", "C3","C4","C5","C6")
    colnames(logp_sparkx) <-c("C1", "C2", "C3","C4","C5","C6")
    rownames(anno_col) <- colnames(logp_ctsv)
    
    png(filename = "text.png",height = 500, width = 600, res = 200)
    p2 <- pheatmap(logp_ctsv[1:1000,],
                   color = colors,
                   breaks = seq(0,max(logp_ctsv[1:1000,]),length.out = length(colors)+1),
                   cluster_cols = F, cluster_rows = F,
                   show_rownames = T, show_colnames = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_cor,
                   annotation_names_col = F,
                   legend_breaks = seq(0,15,5),
                   legend_labels = c("0","5","10","15"),
                   labels_row = labels_row,
                   annotation_legend = F,
                   fontsize_row = 6,
                   width = 3, height = 2.5,
                   filename = paste0("figures/gamma_ctsv_",pat,".png")

    )
    
    p2 <- pheatmap(logp_sparkx[1:1000,],
                   color = colors,
                   cluster_cols = F, cluster_rows = F,
                   show_rownames = T, show_colnames = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_cor,
                   annotation_names_col = F,
                   labels_row = labels_row,
                   fontsize_row = 6,
                   annotation_legend = F,
                   width = 3, height = 2.5,
                   filename = paste0("figures/gamma_sparkx_",pat,".png")
                   
    )
    
    
}
