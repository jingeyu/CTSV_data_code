rm(list = ls())
setwd("~/Desktop/SIM/")
library(SPARK)
library(MASS)
library(DIRECT)


# ----SPARK SE genes----
pattern <- c("linear","gau","cos")
K = 6
tpr_sx2 <- matrix(NA,length(pattern),K)
fp2 <- matrix(NA,length(pattern),K)
for(pat in pattern){
    load(paste0("data/Sim_",pat,".RData"))
    rownames(spot.coor) <- colnames(Y)
    S <- t(spot.coor) - colMeans(spot.coor)
    S <- t(S / apply(S, 1, sd))
    rownames(S) <- colnames(Y)
    G <- nrow(Y)
    n <- ncol(Y)
    seed <- 20210810
    set.seed(seed)
    # sample W hat
    alpha0 <- 100
    W_est <- matrix(NA, K, n)
    for(i in 1:n){
        W_est[,i] <- rDirichlet(1,alpha0 * W[,i])
    }
    indx <- c()
    for(i in 1:n){
        indx[i] <- which.max(W_est[,i]) 
    }
    
    for(k in 1:K){
        indx_k <- which(indx == k)
        subY <- Y[,indx_k]
        subS <- S[indx_k,]
        sparkx.b <- sparkx(subY,subS,numCores=8,option="mixture")
        res1 <- sparkx.b$res_mtest
        write.csv(res1, file = paste0("re/sparkX_sim_",pat,"_ct_",k,".csv"))
        indx_true <- which(gamma_true[,k] == 1)
        indx_se2 <- which(res1[,2] < 0.01)
        tpr_sx2[match(pat,pattern),k] <- length(intersect(indx_se2, indx_true)) / length(indx_true)
        fp2[match(pat,pattern),k] <- length(setdiff(indx_se2,indx_true))
    }
    
}
write.csv(fp2,"re/sparkx_ct_fp_1.csv")
write.csv(tpr_sx2,"re/sparkx_ct_tpr_1.csv")
print(fp2)
print(tpr_sx2)


