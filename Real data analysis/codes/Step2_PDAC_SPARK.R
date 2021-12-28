rm(list = ls())
setwd("")# set your working directory properly
library(SPARK)
st <- 1
load(paste0("data/Real_A_author_st",st,".RData"))
rownames(spot.coor) <- colnames(Y)
spark.st.a <- CreateSPARKObject(counts=Y,
                                location=spot.coor,
                                percentage = 0,
                                min_total_counts = 0)
## total counts for each cell/spot
spark.st.a@lib_size <- apply(spark.st.a@counts, 2, sum)

set.seed(1)
## Estimating Parameter Under Null
t1 <- proc.time()
spark.st.a <- spark.vc(spark.st.a,
                       covariates = NULL,
                       lib_size = spark.st.a@lib_size,
                       num_core = 8,
                       verbose = F)
## Calculating pval
spark.st.a <- spark.test(spark.st.a,
                         check_positive = T,
                         verbose = F)
t2 <- proc.time()
print(t2- t1)

res <- spark.st.a@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
# res <- res[order(res[,2]),]
write.csv(res, file = paste0("re/spark_SEgene_A_ST",st,".csv"))

