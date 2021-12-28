rm(list = ls())
setwd("")# set your working directory properly
library(SPARK)
library(stringr)

st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))
n <- ncol(Y)
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
print(sum(rowSums(W) == 0))

W <- W / rowSums(W)

K <- ncol(W)
print(K)

indx <- apply(W,1,which.max)
print(table(indx))

for(k in 1:K){
    set.seed(1)
    indx_k <- which(indx == k)
    subY <- Y[,indx_k]
    subS <- spot.coor[indx_k,]
    sparkx.b <- sparkx(subY,subS,numCores=8,option="mixture")
    res1 <- sparkx.b$res_mtest
    write.csv(res1, file = paste0("re/sparkX_ZINB_SEgene_A_ST_",st,"ct_",k,".csv"))
    
}