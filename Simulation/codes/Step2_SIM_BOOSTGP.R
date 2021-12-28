rm(list = ls())
setwd("~/Desktop/SIM/")

# ---- Run BOOST-GP ----
source("codes/boost.gp.R")
for(pat in c("linear","gau","cos")){
    load(paste0("data/Sim_",pat,".RData"))
    S <- t(spot.coor) - colMeans(spot.coor)
    S <- t(S / apply(S, 1, sd))
    seed <- 20210810
    set.seed(seed)
    t1 <- proc.time()
    res <- boost.gp(Y = t(Y), loc = S)
    t2 <- proc.time()
    print(t2-t1)
    write.csv(res, file = paste0("re/BOOST_sim_",pat,".csv"))
}
