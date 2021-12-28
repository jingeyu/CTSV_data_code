rm(list = ls())
setwd("~/Desktop/SIM/")
library(SPARK)

# ---- Run SPARK ----
for(pat in c("linear","gau","cos")){
        set.seed(123)
        load(paste0("data/Sim_",pat,".RData"))
        S <- t(spot.coor) - colMeans(spot.coor)
        S <- t(S / apply(S, 1, sd))
        colnames(S) <- c("x","y")
        total_counts <- colSums(Y)
        info <- cbind.data.frame(S, total_counts)
        colnames(info) <- c("x", "y", "UMIS")
        rownames(info) <- colnames(Y)
        
        spark.st <- CreateSPARKObject(counts=Y,
                                      location=info[,1:2],
                                      percentage = 0,
                                      min_total_counts = 0)
        ## total counts for each cell/spot
        spark.st@lib_size <- apply(spark.st@counts, 2, sum)
        t1 <- proc.time()
        spark.st <- spark.vc(spark.st,
                             covariates = NULL,
                             lib_size = spark.st@lib_size,
                             num_core = 8,
                             verbose = F)
        ## Calculating pval
        spark.st <- spark.test(spark.st,
                               check_positive = T,
                               verbose = F)
        t2 <- proc.time()
        print(t2 - t1)
        res <- spark.st@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
        write.csv(res, file = paste0("re/spark_sim_",pat,".csv"))

    }
    

