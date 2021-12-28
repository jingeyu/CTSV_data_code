rm(list = ls())
setwd("")# set your working directory properly
library(pscl)
library(doSNOW)
library(doParallel)
library(stringr)
library(foreach)
library(qvalue)

st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))

#------CTSV------
  G <- nrow(Y)
  cell_type <- colnames(W)
  cell_type <- gsub("[.]","_",cell_type)
  cell_type
  colnames(W) <- cell_type
  indx_c <- which(str_detect(cell_type,"Cancer_clone"))
  indx_m <- which(str_detect(cell_type,"Macrophages"))
  tmp1 <- apply(W[,indx_c], 1, sum)
  tmp2 <- apply(W[,indx_m], 1, sum)
  W <- W[,-c(indx_c,indx_m)]
  W <- cbind(W,tmp1,tmp2)
  colnames(W) <- c(colnames(W)[1:(ncol(W) - 2)],"Cancer clone", "Macrophages")
  wq <- apply(W,2,quantile,0.8)
  indx_out <- which(wq > 0.10)
  W <- W[,indx_out]
  W <- W / rowSums(W)
  K <- ncol(W)
  print("cell type number after filtering")
  print(K)

  print(dim(Y))

  S <- t(spot.coor) - colMeans(spot.coor)
  S <- t(S / apply(S, 1, sd))
  quan <- c(0.2,0.4,0.6,0.8,1.0)
  psi1 <- quantile(abs(S[,1]), quan)
  psi2 <- quantile(abs(S[,2]), quan)

  cl<- makeCluster(8)
  registerDoSNOW(cl)
  pattern <- c("linear","gau1","gau2","cos1","cos2")
  P_VAL <- array(NA, dim = c(G,2*K,5))
  for(pat in pattern){
    if(pat == "gau1"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
    }else if(pat == "gau2"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[3]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[3]^2)
    }else if(pat == "cos1"){
      h1 <- cos(2*pi*S[,1] / psi1[2])
      h2 <- cos(2*pi*S[,2] / psi2[2])
    } else if(pat == "cos2"){
      h1 <- cos(2*pi*S[,1] / psi1[3])
      h2 <- cos(2*pi*S[,2] / psi2[3])
    }else{
      h1 <- S[,1]
      h2 <- S[,2]
    }



    Tmp <- cbind(W * h1, W * h2, W)

    colnames(Tmp) <- 1:ncol(Tmp)

    P_gene <- function(g,Y,Tmp,h1,h2){
      y <- Y[g,]
      K <- ncol(Tmp)/3
      ell <- colSums(Y) / median(colSums(Y))
      err <- try(fm_zinb0 <- zeroinfl(y ~ -1+offset(log(ell))+Tmp|1,
                                      dist = "negbin",link = "probit",
                                      control = zeroinfl.control(method = "CG"
                                      )), silent = TRUE)
      if(class(err) == 'try-error') {
        p_val <- rep(-1,2*K)
      } else{
        p_val <- coef(summary(fm_zinb0))$count[,4]
        nind <- 2*(length(p_val) - 1)/3
        p_val <- p_val[1:nind]

      }
      return(p_val)

    }
    
    t1 <- proc.time()
    P_VAL[,,match(pat,pattern)] = foreach(g=1:G,
                                              .combine=rbind,
                                              .packages=c("pscl")) %dopar% P_gene(g,Y,Tmp,h1,h2)
    
    rownames(P_VAL[,,match(pat,pattern)]) <- rownames(Y)
    
    t2 <- proc.time()
    print(t2-t1)
  }
stopCluster(cl)

# Cauchy combination rule
P_VAL[which(is.na(P_VAL))] <- 1
P_VAL[P_VAL == -1] <- 1
P_VAL <- tan((0.5 - P_VAL)*pi)
T_cau0 <- apply(P_VAL, c(1,2), mean)
P_val <- 1-pcauchy(T_cau0)
# convert q-values into q-values
Q_val <- matrix(qvalue(c(P_val))$qvalue,G,2*K)
rownames(Q_val) <- colnames(Y)
    
write.csv(P_val, file = "re/CTSV_A1_pval.csv")
write.csv(Q_val, file = "re/CTSV_A1_qval.csv")

