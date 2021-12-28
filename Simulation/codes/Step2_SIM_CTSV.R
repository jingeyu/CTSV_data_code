rm(list = ls())
setwd("~/Desktop/SIM/")
library(DIRECT)
library(pscl)
library(doSNOW)
library(qvalue)

###############################################################
########################### CTSV ##############################
###############################################################
#--------
load("data/spot_filter_ct4.RData")
# normalize coordinates
S <- t(spot.coor) - colMeans(spot.coor)
S <- t(S / apply(S, 1, sd))

#spot number
n <- nrow(spot.coor)

#gene number
G <- 10000

# cell type number
K <- 6

quan <- c(0.2,0.4,0.6,0.8,1.0)
psi1 <- quantile(abs(S[,1]), quan)
psi2 <- quantile(abs(S[,2]), quan)

#---- Run CTSV for three spatial patterns -----
# set cores for paralleling

for(pat in c("linear","gau","cos")){
    for(fit_pat in c("gau2","cos2","gau1","cos1","linear")){
        load(paste0("data/Sim_",pat,".RData"))
        
        seed <- 20210810
        set.seed(seed)
        # sample W hat
        alpha0 <- 100
        W_est <- matrix(NA, K, n)
        for(i in 1:n){
            W_est[,i] <- rDirichlet(1,alpha0 * W[,i])
        }
        
        if(fit_pat == "gau1"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
        }else if(fit_pat == "gau2"){
            h1 <- exp(-S[,1]^2 / 2 / psi1[3]^2)
            h2 <- exp(-S[,2]^2 / 2 / psi2[3]^2)
        }else if(fit_pat == "cos1"){
            h1 <- cos(2*pi*S[,1] / psi1[2])
            h2 <- cos(2*pi*S[,2] / psi2[2])
        } else if(fit_pat == "cos2"){
            h1 <- cos(2*pi*S[,1] / psi1[3])
            h2 <- cos(2*pi*S[,2] / psi2[3])
        }else{
            h1 <- S[,1]
            h2 <- S[,2]
        }
        
        Tmp <- cbind(t(W_est) * h1, t(W_est) * h2, t(W_est))
        colnames(Tmp) <- 1:ncol(Tmp)
        G <- nrow(Y)
        
        P_gene <- function(g,Y,Tmp,h1,h2){
            y <- Y[g,]
            ell <- colSums(Y) / median(colSums(Y))
            fm_zinb0 <- zeroinfl(y ~ -1+offset(log(ell))+Tmp|1,
                                 dist = "negbin",link = "probit",
                                 control = zeroinfl.control(method = "CG"
                                 ))
            p_val <- coef(summary(fm_zinb0))$count[,4]
            nind <- 2*(length(p_val) - 1)/3
            p_val <- p_val[1:nind]
            write.table(g, file = paste0("iteration_",pat,"_",fit_pat,".txt"),col.names = F,append = T)
            return(p_val)
        }
        
        cl<- makeCluster(8)
        registerDoSNOW(cl)
        t1 <- proc.time()
        P_val = foreach(g=1:G,
                        .combine=rbind,
                        .packages=c("pscl")) %dopar% P_gene(g,Y,Tmp,h1,h2)
        
        t2 <- proc.time()
        print(t2-t1)
        stopCluster(cl)
        rownames(P_val) <- rownames(Y)
        write.csv(P_val, file = paste0("re/SE_CTSV_CG_",pat,"_",fit_pat,"_p.csv"))
    }
    
}


# transform p-values into combined p-values
for(pat in c("linear", "gau", "cos")){
    P_ZINB1 <- read.csv(file = paste0("re/SE_CTSV_CG_",pat,"_","linear","_p.csv"))
    rownames(P_ZINB1) <- P_ZINB1[,1]
    P_ZINB1 <- as.matrix(P_ZINB1[,-1])
    P_ZINB1[which(is.na(P_ZINB1))] <- 1
    
    G <- nrow(P_ZINB1)
    K <- ncol(P_ZINB1)/2
    
    P_val_ZINB_cos1 <- read.csv(file = paste0("re/SE_CTSV_CG_",pat,"_","cos1","_p.csv"))
    rownames(P_val_ZINB_cos1) <- P_val_ZINB_cos1[,1]
    P_val_ZINB_cos1 <- as.matrix(P_val_ZINB_cos1[,-1])
    P_val_ZINB_cos2 <- read.csv(file = paste0("re/SE_CTSV_CG_",pat,"_","cos2","_p.csv"))
    rownames(P_val_ZINB_cos2) <- P_val_ZINB_cos2[,1]
    P_val_ZINB_cos2 <- as.matrix(P_val_ZINB_cos2[,-1])
    P_val_ZINB_gau1 <- read.csv(file = paste0("re/SE_CTSV_CG_",pat,"_","gau1","_p.csv"))
    rownames(P_val_ZINB_gau1) <- P_val_ZINB_gau1[,1]
    P_val_ZINB_gau1 <- as.matrix(P_val_ZINB_gau1[,-1])
    P_val_ZINB_gau2 <- read.csv(file = paste0("re/SE_CTSV_CG_",pat,"_","gau2","_p.csv"))
    rownames(P_val_ZINB_gau2) <- P_val_ZINB_gau2[,1]
    P_val_ZINB_gau2 <- as.matrix(P_val_ZINB_gau2[,-1])

    P_val_ZINB_cos1[which(is.na(P_val_ZINB_cos1))] <- 1
    P_val_ZINB_cos2[which(is.na(P_val_ZINB_cos2))] <- 1
    P_val_ZINB_gau1[which(is.na(P_val_ZINB_gau1))] <- 1
    P_val_ZINB_gau2[which(is.na(P_val_ZINB_gau2))] <- 1
    
    #---------Cauchy combination rule---------
    
    T_cau0 <- (tan((0.5-P_ZINB1) * pi) + tan((0.5-P_val_ZINB_gau1) * pi) + 
                   tan((0.5-P_val_ZINB_cos1) * pi)+
                   tan((0.5-P_val_ZINB_gau2) * pi) + 
                   tan((0.5-P_val_ZINB_cos2) * pi) ) / 5
    P_val_ZINB <- 1-pcauchy(T_cau0)
    Q_ZINB1 <- matrix(qvalue(c(P_val_ZINB))$qvalue,G,2*K)
    
    P_val_ZINB1<- matrix(NA,G,K)
    Q_val_ZINB1<- matrix(NA,G,K)
    for(k in 1:K){
        P_val_ZINB1[,k] <- apply(P_val_ZINB[,c(k,k+K)],1,min)
        Q_val_ZINB1[,k] <- apply(Q_ZINB1[,c(k,k+K)],1,min)
    }
    rownames(P_val_ZINB) <- rownames(P_ZINB1)
    rownames(Q_ZINB1) <- rownames(P_ZINB1)
    
    write.csv(P_val_ZINB, file = paste0("re/CTSV_sim_",pat,"_pval.csv"))
    write.csv(Q_ZINB1, file = paste0("re/CTSV_sim_",pat,"_qval.csv"))
    
}





