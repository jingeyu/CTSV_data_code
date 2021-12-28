rm(list = ls())
setwd("~/Desktop/SIM")
library(MASS)
library(DIRECT)

###############################################################
###################### Data Generation ########################
###############################################################
#---- load spatial coordinates ----
load("data/spot_filter_ct4.RData")
# normalize coordinates
S <- t(spot.coor) - colMeans(spot.coor)
S <- t(S / apply(S, 1, sd))

#spot number
n <- nrow(spot.coor)

#gene number
G <- 10000

# cell type number
K <- 6

# number of DE gene
DE_num <- 300

# size parameter in NB distribution
size <- 100

# zero inflation rate
pai0 <- 0.60

# ## ---generate data-----

for(pat in c("linear","gau","cos")){
        seed <- 20210509
        set.seed(seed)
        NDE_scrna <- rnorm(G, mean=2, sd=0.2)
        scrna_1 <- NDE_scrna
        scrna_2 <- NDE_scrna
        scrna_3 <- NDE_scrna
        scrna_4 <- NDE_scrna
        scrna_5 <- NDE_scrna
        scrna_6 <- NDE_scrna
        scrna_2[sample(1:G,DE_num,replace = F)] <- rnorm(DE_num, mean=3, sd=0.2)
        scrna_3[sample(1:G,DE_num,replace = F)] <- rnorm(DE_num, mean=2, sd=0.2)
        scrna_4[sample(1:G,DE_num,replace = F)] <- rnorm(DE_num, mean=4, sd=0.2)
        scrna_5[sample(1:G,DE_num,replace = F)] <- rnorm(DE_num, mean=3, sd=0.2)
        scrna_6[sample(1:G,DE_num,replace = F)] <- rnorm(DE_num, mean=4, sd=0.2)
        
        eta <- cbind(scrna_1,scrna_2,scrna_3,scrna_4,scrna_5,scrna_6)
        
        gamma_true <- matrix(0, G, K)
        gamma_true[101:300,1] <- 1
        gamma_true[201:400,2] <- 1
        gamma_true[301:500,3] <- 1
        gamma_true[401:600,4] <- 1
        gamma_true[501:700,5] <- 1
        gamma_true[601:800,6] <- 1
        
        beta1 <- matrix(0, G, K)
        beta2 <- matrix(0, G, K)
        
        alpha1 <- c(1,1,1,1,1,1)
        alpha2 <- c(1,3,5,7,9,11)
        alpha3 <- c(16,14,12,10,8,6)
        alpha4 <- c(1,4,4,4,4,1)
        alpha <- cbind(alpha1, alpha2, alpha3, alpha4)
        
        W <- matrix(NA, K, n)
        for(i in 1:n){
            W[, i] <- rDirichlet(1, alpha[, spot.region[i]])
        }
        
        if(pat == "linear"){
            h1 <- S[,1]
            h2 <- S[,2]
            beta1[gamma_true == 1] <- 1.8
            beta2[gamma_true == 1] <- 0.8
        }else if(pat == "gau"){
            h1 <- exp(-S[,1]^2 / 2)
            h2 <- exp(-S[,2]^2 / 2)
            beta1[gamma_true == 1] <- 3
            beta2[gamma_true == 1] <- 1
        }else{
            h1 <- cos(2*pi*S[,1])
            h2 <- cos(2*pi*S[,2])
            beta1[gamma_true == 1] <- 2.5
            beta2[gamma_true == 1] <- 1
        }
        
        
        log_lambda <- eta %*% W + beta1 %*% t(t(W) * h1) + beta2 %*% t(t(W) * h2)
        Y <- matrix(rnbinom(G*n,size = size, mu = exp(c(log_lambda))), G, n)
        set.seed(5)
        r_unif <- matrix(runif(G*n),G,n)
        Y[r_unif <= pai0] <- 0
        # print(summary(c(Y)))
        
        rownames(Y) <- paste0("Gene",1:G)
        colnames(Y) <- paste0("Spot",1:n)
        
        save(list = c("Y", "eta", "h1", "h2", "beta1", "beta2",
                      "W", "spot.coor", "gamma_true","log_lambda"),
             file = paste0("data/Sim_",pat,".RData"))
        
}


for(pat in c("linear","cos","gau")){
                load(paste0("data/Sim_",pat,".RData"))
                write.csv(t(Y),file = paste0("csvs/Sim_",pat,".csv"))
                total_counts <- colSums(Y)
                st_a_info <- cbind(S,total_counts)
                rownames(st_a_info) <- colnames(Y)
                colnames(st_a_info) <- c("x","y","total_counts")
                write.csv(st_a_info,file = paste0("csvs/Sim_info_",pat,".csv"))

}

