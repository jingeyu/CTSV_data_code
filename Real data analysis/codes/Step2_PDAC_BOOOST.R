rm(list = ls())
setwd("")# set your working directory properly
library(MASS)
# library(pROC)
###############################################################
###################### BOOST-GP ########################
###############################################################
source("codes/boost.gp.R")

st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))
    
set.seed(1)
t1 <- proc.time()
res <- boost.gp(Y = t(Y), loc = as.data.frame(spot.coor))
t2 <- proc.time()
print(t2-t1)
write.csv(res, file = paste0("re/BOOST_res_A_",st,".csv"))

