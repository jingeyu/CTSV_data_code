rm(list = ls())
setwd("~/Desktop/SIM")
library(dplyr)
###############################################################
##################### Teble 1 and table2 ######################
###############################################################
load("data/Sim_linear.RData")
G <- nrow(Y)
K <- nrow(W)
pattern <- c("linear","gau","cos")
true_se <- rep(0,G)
true_se[which(rowSums(gamma_true) > 0)] <- 1
#-------TPR + FP -------
# We didd not include trendsceekad somde since they did not identify any SV gene
alpha <- 0.01
pai0 <- 0.6
#----1------
TPR_BOOST <- matrix(NA, length(pattern),length(pai0))
TPR_SPARK <- matrix(NA, length(pattern),length(pai0))
TPR_Spatial <- matrix(NA, length(pattern),length(pai0))
TPR_ZINB_T <- matrix(NA, length(pattern),length(pai0))
TPR_SPARKX_T <- matrix(NA, length(pattern),length(pai0))

FP_BOOST <- matrix(NA, length(pattern),length(pai0))
FP_SPARK <- matrix(NA, length(pattern),length(pai0))
FP_Spatial <- matrix(NA, length(pattern),length(pai0))
FP_ZINB_T <- matrix(NA, length(pattern),length(pai0))
FP_SPARKX_T <- matrix(NA, length(pattern),length(pai0))
TPR_CT <- NULL
FP_CT <- NULL
num <- 1
#-----2----
for(pat in pattern){
    
        Q_ZINB1 <- read.csv(paste0("re/CTSV_sim_",pat,"_qval.csv"))
        rownames(Q_ZINB1) <- Q_ZINB1[,1]
        Q_ZINB1 <- Q_ZINB1[,-1]
        P_val_ZINB1<- matrix(NA,G,K)
        Q_val_ZINB1 <- matrix(NA,G,K)
        for(k in 1:K){
          
            Q_val_ZINB1[,k] <- apply(Q_ZINB1[,c(k,k+K)],1,min)
        }
        rownames(Q_val_ZINB1) <- rownames(Q_ZINB1)
        q_val_zinb1 <- apply(Q_val_ZINB1,1,min)
        zinb_gene <- rownames(Q_val_ZINB1)
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
        
        SE_spark <- spark_se[,1][spark_se[,3] < alpha]
        SE_sparkx <- sparkx_se[,1][sparkx_se[,3] < alpha]
        SE_spatial <- spatial_gene[spatial_se[,3] < alpha]
        SE_ZINB2 <- zinb_gene[q_val_zinb1 < alpha]
        SE_boost <- rownames(boost_se)[boost_ppi > 0.5]
        SE_true <- rownames(Y)[true_se == 1]
        TPR_SPARK[match(pat,pattern),num] <- length(intersect(SE_spark, SE_true))/ sum(true_se)
        TPR_SPARKX_T[match(pat,pattern),num] <- length(intersect(SE_sparkx, SE_true))/ sum(true_se)
        TPR_Spatial[match(pat,pattern),num] <- length(intersect(SE_spatial, SE_true))/ sum(true_se)
        TPR_ZINB_T[match(pat,pattern),num] <- length(intersect(SE_ZINB2, SE_true))/ sum(true_se)
        TPR_BOOST[match(pat,pattern),num] <- length(intersect(SE_boost, SE_true))/ sum(true_se)
        
        FP_SPARK[match(pat,pattern),num] <- length(setdiff(SE_spark, SE_true))
        FP_SPARKX_T[match(pat,pattern),num] <- length(setdiff(SE_sparkx, SE_true))
        FP_Spatial[match(pat,pattern),num] <- length(setdiff(SE_spatial, SE_true))
        FP_ZINB_T[match(pat,pattern),num] <- length(setdiff(SE_ZINB2, SE_true))
        FP_BOOST[match(pat,pattern),num] <- length(setdiff(SE_boost, SE_true))

        
        FP_zinb <- c()
        TPR_zinb <- c()
        TPR_spk <- c()
        FP_spk <- c()
        
        for(k in 1:K){
            spk <- read.csv(file = paste0("re/sparkX_sim_",pat,"_ct_",k,".csv"))
            SE_true_k <- rownames(Y)[gamma_true[,k]==1]
            positive_total <- length(SE_true_k)
            spk_se <- spk$X[spk$adjustedPval < alpha]
            ZINB_se <- zinb_gene[Q_val_ZINB1[,k] < alpha]
            TPR_zinb[k] <-  length(intersect(SE_true_k, ZINB_se)) / positive_total
            FP_zinb[k] <- length(setdiff(ZINB_se, SE_true_k))
            TPR_spk[k] <-  length(intersect(SE_true_k, spk_se)) / positive_total
            FP_spk[k] <- length(setdiff(spk_se, SE_true_k))
            
        }
        
        TPR_ct <- cbind(TPR_zinb,TPR_spk)
        TPR_CT <- cbind(TPR_CT,TPR_ct)
        FP_ct <- cbind(FP_zinb,FP_spk)
        FP_CT <- cbind(FP_CT,FP_ct)
    }



#----3-----

rownames(FP_CT) <- paste0("cell-type ",1:K)
rownames(TPR_CT) <- paste0("cell-type ",1:K)
colnames(FP_CT) <- rep(c("CTSV","SPARKX"),length(pattern))
colnames(TPR_CT) <- rep(c("CTSV","SPARKX"),length(pattern))

write.csv(TPR_CT, file = paste0("result/TPR_sim_ct_",100*alpha,".csv"))
write.csv(FP_CT, file = paste0("result/FR_sim_ct_",100*alpha,".csv"))


TPR <- rbind(TPR_SPARK,TPR_SPARKX_T,TPR_Spatial,TPR_ZINB_T,TPR_BOOST)
FP <- rbind(FP_SPARK,FP_SPARKX_T,FP_Spatial,FP_ZINB_T,FP_BOOST)

com_me <- c("SPARK","SPARKX","SpatislDE","CTSV","BOOST")
result <- cbind(TPR,FP)
rownames(result) <- lapply(com_me,function(x) paste(x,combn(pattern,1),sep = "_")) %>% unlist
colnames(result) <- c("TPR","FP")
write.csv(result, file = paste0("result/com_methods",100*alpha,".csv"))


