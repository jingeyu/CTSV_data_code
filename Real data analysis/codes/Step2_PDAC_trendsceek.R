rm(list = ls())
setwd("")# set your working directory properly
library(trendsceek)
source("codes/trendsceek.R")
###############################################################
###################### trendsceek ########################
###############################################################
st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))

set.seed(1)

pp <- pos2pp(spot.coor)
##Set marks as the logged normalized gene expression
min.count = 1
counts_norm = deseq_norm(Y, min.count)
log.fcn = log10
pp = set_marks(pp, counts_norm, log.fcn = log.fcn)


#=======Run trendsceek========
##set parameters
nrand = 10
ncores = 1

t1 <- proc.time()
trendstat_list = trendsceek_test(pp, nrand, ncores)  
t2 <- proc.time()

##p.BH
supstats_list = trendstat_list[['supstats']]
G <- nrow(supstats_list[[1]])
P_val <- matrix(NA,G,4)
Genes <- list()
for(i in 1:4){
    P_val[,i] <- supstats_list[[i]][["p.bh"]]
    Genes[[i]] <- supstats_list[[i]][["gene"]]
}
for(i in 1:4){
    tmp <- match(rownames(Y),Genes[[i]])
    P_val[,i] <- P_val[,i][tmp]
}
rownames(P_val) <- rownames(Y)
colnames(P_val) <- names(supstats_list)
print(sum(P_val) < 0.05)
write.csv(P_val, file = paste0("re/trend_A_",st,".csv"))
