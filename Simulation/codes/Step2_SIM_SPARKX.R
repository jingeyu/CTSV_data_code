rm(list = ls())
setwd("~/Desktop/SIM/")
library(SPARK)

# ---- Run SPARKX ----
pattern <- c("gau","cos","linear")
for(pat in pattern){
    load(paste0("data/Sim_",pat,".RData"))
    S <- t(spot.coor) - colMeans(spot.coor)
    S <- t(S / apply(S, 1, sd))
    rownames(S) <- colnames(Y)
    set.seed(123)
    t1 <- proc.time()
    sparkx.b <- sparkx(Y,S,numCores=8,option="mixture")
    t2 <- proc.time()
    print(t2 - t1)
    res1 <- sparkx.b$res_mtest
    write.csv(res1, file = paste0("re/sparkx_sim_",pat,".csv"))
   
}
