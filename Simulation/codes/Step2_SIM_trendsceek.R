rm(list = ls())
setwd("~/Desktop/SIM/")
library(trendsceek)

source("codes/trendsceek.R")

for(pat in c("linear","gau","cos")){
    load(paste0("data/Sim_",pat,".RData"))
    S <- t(spot.coor) - colMeans(spot.coor)
    t1 <- proc.time()
    S <- t(S / apply(S, 1, sd))
    ##Convert cell positions to point pattern
    pp = pos2pp(S)
    ##Set marks as the logged normalized gene expression
    min.count = 1
    counts_norm = deseq_norm(Y, min.count)
    log.fcn = log10
    pp = set_marks(pp, counts_norm, log.fcn = log.fcn)

    #=======Run trendsceek========
    ##set parameters
    nrand = 10
    ncores = 8

    set.seed(123)
    trendstat_list = trendsceek_test(pp, nrand, ncores)  
    t2 <- proc.time()
    print(t2-t1)
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
    write.csv(P_val, file = paste0("re/trend_sim_",pat,".csv"))
    
    
}
