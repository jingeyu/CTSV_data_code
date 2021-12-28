rm(list = ls())
setwd("")# set your working directory properly
library(SPARK)
st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))
rownames(spot.coor) <- colnames(Y)
set.seed(1)
t1 <- proc.time()
sparkx.b <- sparkx(Y,spot.coor,numCores=6,option="mixture")
write.csv(sparkx.b$res_mtest, file = paste0("re/sparkX_SEgene_A_ST",st,"_all.csv"))
t2 <- proc.time()
print(t2-t1)
