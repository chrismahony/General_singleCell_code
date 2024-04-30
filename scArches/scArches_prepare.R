
#prepare amp2 data
load("/rds/projects/c/croftap-celldive01/amp2/amp2_full.RData")
Idents(amp2_new_f)<- 'cluster_name'
fibs_amp2 <- subset(amp2_new_f, idents=levels(amp2_new_f)[grep("F-", levels(amp2_new_f))])
rm(list=ls()[! ls() %in% c("fibs_amp2")])
gc()
setwd("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches")

fibs_amp2 <- fibs_amp2 %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims = 1:30)

library(SeuratDisk)
SaveH5Seurat(fibs_amp2, filename = "/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_amp2_seurat.h5Seurat")
Convert("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_amp2_seurat.h5Seurat", dest = "h5ad")


#prepare MAPJAG data
load("/rds/projects/c/croftap-mapjagdata/MAPJAGv2/Chris/CM_analysis.RData")
DimPlot(stroma_clean_h)
rm(list=ls()[! ls() %in% c("fibs_amp2", "stroma_clean_h")])
gc()


stroma_clean_h_meta <- stroma_clean_h@meta.data

stroma_clean_h$clusters <- as.character(stroma_clean_h$clusters)

stroma_clean_h$global_new <- as.character(stroma_clean_h$global_new)
stroma_clean_h$TYPE <- as.character(stroma_clean_h$TYPE)
stroma_clean_h$condition <- as.character(stroma_clean_h$condition)
stroma_clean_h$sample <- as.character(stroma_clean_h$sample)
stroma_clean_h$batch <- as.character(stroma_clean_h$batch)


Idents(stroma_clean_h) <- 'clusters'
fibs_mapjag <- subset(stroma_clean_h, idents=levels(stroma_clean_h)[c(1,3,4,5,7,8,10)])

fibs_mapjag <- fibs_mapjag %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)



stroma_clean_h_meta <- stroma_clean_h@meta.data


SaveH5Seurat(stroma_clean_h, filename = "/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/stroma_clean_h_seurat.h5Seurat")
Convert("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/stroma_clean_h_seurat.h5Seurat", dest = "h5ad")

SaveH5Seurat(fibs_mapjag, filename = "/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_mapjag_seurat.h5Seurat", overwrite = T)
Convert("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/fibs_mapjag_seurat.h5Seurat", dest = "h5ad", overwrite = T)

