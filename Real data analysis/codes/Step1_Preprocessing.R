rm(list = ls())
setwd("")
# set your working directory properly

#=== Preprocess =====
pdac_A <- readRDS(file = "processed_data/PDAC-A_itai_joint.RDS")
set.seed(1)
pdac_A <- Seurat::SCTransform(object = pdac_A)
pdac_A <- Seurat::RunPCA(pdac_A, verbose = FALSE)
Seurat::ElbowPlot(pdac_A, ndims = 50)

pdac_A <- Seurat::FindNeighbors(pdac_A,
                                dims = 1:40)
pdac_A <- Seurat::FindClusters(pdac_A,
                               verbose = FALSE, 
                               resolution = c(1, 2, 3, 4, 5))
pdac_A <- Seurat::RunUMAP(pdac_A,
                          dims = 1:40)

pdac_a_annot <- data.frame(
    annotation = sort(unique(as.character(pdac_A$annotation))),
    new_name = c("Acinar cells", "Cancer clone (TM4SF1)", "Cancer clone (S100A4)", 
                 "High/Hypoxic ductal cells (APOL1)", "Centroacinar ductal cells",
                 "Antigen-presenting ductal cells", 
                 "Terminal ductal cells", "Endocrine cells", "Endothelial cells", 
                 "Fibroblasts", "Macrophages M2", "Macrophages M1", "Mast cells",
                 "mDCs", "mDCs", "Monocytes", "pDCs", "RBCs", "T cells & NK cells",
                 "Tuft cells"))

new_annot <- data.frame(pdac_A@meta.data) %>% 
    left_join(pdac_a_annot, by = "annotation") %>% 
    pull(new_name) %>% 
    as.character()

pdac_A@meta.data[["annotation"]] <- new_annot
Idents(pdac_A) <- pdac_A$annotation

saveRDS(object = pdac_A,
        file = "processed_data/PDAC-A_itai_processed.RDS")


#=====Deconvolution=======
tech <- "indrop"
tissue <- "pdac_itai_2"
dwn_smplng <- "both"
org <- "hs"

clust_vr <- "annotation"
cl_n <- 100
method <- "nsNMF"
transf <- "uv"
hvg <- 3000
FC <- 1
pct1 <- 0.9

options(stringsAsFactors = FALSE)


indrop_pdac_a <- readRDS(file = "processed_data/PDAC-A_itai_processed.RDS")
dim(indrop_pdac_a)

# Remove RBCs
indrop_pdac_a <- indrop_pdac_a[, indrop_pdac_a$annotation != "RBCs"]
dim(indrop_pdac_a)

length(table(Idents(indrop_pdac_a)))

# Spatial deconvolution

# Find markers in each cluster
Seurat::Idents(object = indrop_pdac_a) <- indrop_pdac_a@meta.data[, clust_vr]
cluster_markers_a <- Seurat::FindAllMarkers(object = indrop_pdac_a,
                                            verbose = TRUE,
                                            only.pos = TRUE,
                                            assay = "SCT",
                                            slot = "data")

cluster_markers_filt_a <- cluster_markers_a %>%
    filter(avg_log2FC > 1 & pct.1 > 0.75)

cluster_markers_filt_a$cluster <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".", 
                                       x = cluster_markers_filt_a$cluster, 
                                       perl = TRUE)
saveRDS(cluster_markers_filt_a, file = "processed_data/Markersall_padc_A.RDS")

# Spatial deconvolution
set.seed(1)
st_a_list <- readRDS("processed_data/PDAC-A_ST_list.RDS")

i <- 1
st_a_se <- st_a_list[[i]]

decon_mtrx_a <- SPOTlight::spotlight_deconvolution(se_sc = indrop_pdac_a,
                                                   counts_spatial = st_a_se@assays$RNA@counts,
                                                   cluster_markers = cluster_markers_filt_a,
                                                   cl_n = cl_n,
                                                   hvg = hvg,
                                                   ntop = NULL,
                                                   transf = transf,
                                                   clust_vr = clust_vr,
                                                   method = method,
                                                   min_cont = 0.01)


saveRDS(object = decon_mtrx_a,
        file = paste0("data/STA_ST",i,"_decon.RDS"))


  sp1 <- colnames(st_a_se)
  sp1 <- str_split(sp1,"_")
  index_get <- function(x){
    sp2 <- str_split(x[2],"x")[[1]]
    return(as.numeric(sp2))
  }
  spot.coor <- t(sapply(sp1,index_get))
  W <- readRDS(paste0("data/STA_ST",i,"_decon.RDS"))
  W <- W[[2]][,-ncol(W[[2]])]
  spot.coor <- t(sapply(sp1,index_get))
  st_a_se <- CreateSeuratObject(counts = st_a_se[["RNA"]]@counts,
                                min.cells = 20,
                                min.features = 1)
  Y <- as.matrix(st_a_se[["RNA"]]@counts)
  K <- ncol(W)
  save(list = c("W", "K", "Y","spot.coor"),
       file = paste0("data/Real_A_author_st",i,".RData"))


#----- transform to csv files for the sake of SpatialDE and SOMDE analysis----
load(paste0("Real_A_author_st",i,".RData"))
write.csv(t(Y),file = paste0("data/PDAC_A_counts_",i,".csv"))
total_counts <- colSums(Y)
st_a_info <- cbind(spot.coor,total_counts)
rownames(st_a_info) <- colnames(Y)
colnames(st_a_info) <- c("x","y","total_counts")
write.csv(st_a_info,file = paste0("data/PDAC_A_info_",i,".csv"))


