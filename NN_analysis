
#all is my seurat obj from qpath segmentation annotated using makrer genes with all sample in the experiment

all_obj_split <- SplitObject(all, split.by = "orig.ident")

data_x_y_niches <- data_x_y
library(dbscan)
eps <- 50
nn_list <- list()
nn_df_list <- list()
cluster_ids <- list()
nn_count <- list()
nn_mat <- list()
k_means_res <- list()
k_means_id <- list()
k_means_df <- list()
nn_objs <- list()

for (i in 1:length(data_x_y_niches)){
  nn_list[[i]] <- frNN(x= data_x_y_niches[[i]] %>% as.matrix(), eps = eps)
  nn_df_list[[i]]<- nn_list[[i]]$id %>%
  stack()
cluster_ids[[i]] <- all_obj_split[[i]]$named_niches_final %>% unname()
nn_df_list[[i]]$cluster_id<- cluster_ids[[i]][nn_df_list[[i]]$values]
nn_df_list[[i]]$cluster_id<- factor(nn_df_list[[i]]$cluster_id)
nn_count[[i]]<- nn_df_list[[i]] %>%
  group_by(ind) %>%
  dplyr::count(cluster_id, .drop=F)
nn_count[[i]]<- nn_count[[i]] %>%
  tidyr::pivot_wider(names_from = cluster_id, values_from = n)
nn_mat[[i]]<- nn_count[[i]][,-1] %>% as.matrix()
rownames(nn_mat[[i]])<- nn_count[[i]]$ind
k_means_res[[i]]<- kmeans(nn_mat[[i]], centers = 6)
  k_means_id[[i]]<- k_means_res[[i]]$cluster %>%
  tibble::enframe(name = "cell_id", value = "kmeans_cluster")
  k_means_df[[i]]<- as.data.frame(k_means_id[[i]])
rownames(k_means_df[[i]])<- k_means_id[[i]]$cell_id
  nn_objs[[i]] <- CreateSeuratObject(counts = t(nn_mat[[i]]),  min.features = 1)
  DefaultAssay(nn_objs[[i]]) <- 'RNA'
  nn_objs[[i]] <- ScaleData(nn_objs[[i]])
  nn_objs[[i]] <- RunPCA(nn_objs[[i]], npcs = 10, features = rownames(nn_objs[[i]]))
  nn_objs[[i]] <- FindNeighbors(nn_objs[[i]], reduction = "pca", dims = 1:5)
  nn_objs[[i]] <- RunUMAP(nn_objs[[i]], dims = 1:5)
  nn_objs[[i]] <- FindClusters(nn_objs[[i]], resolution = c(0.02)) 
}

names(nn_objs) <- names(data_x_y_niches)
names_orig <- names(data_x_y_niches)

for (i in 1:length(nn_objs)){
nn_objs[[i]]$orig.ident <- names[[i]] 
}

library(harmony)
all_niches_merged <- merge(x=nn_objs[[1]], y=nn_objs[-1])
all_niches_merged <-ScaleData(all_niches_merged)
all_niches_merged <-FindVariableFeatures(all_niches_merged)
all_niches_merged<-RunPCA(all_niches_merged)
all_niches_merged<-RunUMAP(all_niches_merged, dims = 1:5)
DimPlot(all_niches_merged, group.by="orig.ident", raster=FALSE)
all_niches_merged <- RunHarmony.Seurat_CM(all_niches_merged, group.by.vars = "orig.ident")
all_niches_merged <- RunUMAP(all_niches_merged, reduction="harmony", dims=1:5)
DimPlot(all_niches_merged, group.by = "orig.ident", raster = F)

all_niches_merged <- FindNeighbors(all_niches_merged, reduction="harmony", dims = 1:5)
all_niches_merged <- FindClusters(all_niches_merged, resolution = c(0.02, 0.01, 0.05))
all_niches_merged <- FindClusters(all_niches_merged, resolution = c(0.03)) #try different resolutions, ideally pick one that works for all smaples (this can behard)

DimPlot(all_niches_merged, raster=F, group.by = "RNA_snn_res.0.01") #mine looked really weird
all_niches_merged_slit <- SplitObject(all_niches_merged, split.by = "orig.ident")

library(data.table)
ggplot_niches <- list()
for (i in 1:length(all_niches_merged_slit)){
data_x_y_niches[[i]] <- data_x_y_niches[[i]][rownames(data_x_y_niches[[i]]) %in% colnames(nn_objs[[i]]),]
data_x_y_niches[[i]]$niches <- all_niches_merged_slit[[i]]$RNA_snn_res.0.05
data_x_y_niches[[i]]$cell_id <- colnames(all_niches_merged_slit[[i]])
ggplot_niches[[i]] <- ggplot(data_x_y_niches[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.1)+theme_classic()+ggtitle("niches_spatial") 
print(ggplot_niches[[i]])  
}
table(all_niches_merged$RNA_snn_res.0.03)


#remove clusters with <1000 cells (I  had lost of clusters with about 100-200 cells)
Idents(all_niches_merged) <- 'RNA_snn_res.0.03'
all_niches_merged_fil <- subset(all_niches_merged, idents=c("0", "2", "1", "3", "4", "5"))
all_niches_merged_slit <- SplitObject(all_niches_merged_fil, split.by = "orig.ident")
library(data.table)
data_x_y_niches_f <- data_x_y_niches
ggplot_niches <- list()
smaple_id <- c(1:12)
for (i in 1:length(all_niches_merged_slit)){
rownames(data_x_y_niches_f[[i]]) <- paste(rownames(data_x_y_niches_f[[i]]), smaple_id[[i]], sep="_")
data_x_y_niches_f[[i]] <- data_x_y_niches_f[[i]][rownames(data_x_y_niches_f[[i]]) %in% colnames(all_niches_merged_slit[[i]]),]
data_x_y_niches_f[[i]]$niches <- all_niches_merged_slit[[i]]$RNA_snn_res.0.03
data_x_y_niches_f[[i]]$cell_id <- colnames(all_niches_merged_slit[[i]])
ggplot_niches[[i]] <- ggplot(data_x_y_niches_f[[i]], aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=niches)) +
    geom_point(size=0.1)+theme_classic()+ggtitle("niches_spatial") 
print(ggplot_niches[[i]])  
}

cell_fun = function(j, i, x, y, width, height, fill) {
                grid::grid.rect(x = x, y = y, width = width *0.99, 
                                height = height *0.99,
                                gp = grid::gpar(col = "grey", 
                                                fill = fill, lty = 1, lwd = 0.5))
}

col_fun=circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


avg_abun<- AverageExpression(
  all_niches_merged_fil,
  assays = NULL,
  features = rownames(all_niches_merged_fil),
  return.seurat = FALSE,
  group.by = "RNA_snn_res.0.03")

Heatmap(t(scale(t(avg_abun$RNA))),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = cell_fun,
        col = col_fun,
        column_names_rot = 45)


library(data.table)
data_x_y_niches_all <- rbindlist(data_x_y_niches_f)
all_niches_f <- all[,colnames(all) %in% data_x_y_niches_all$cell_id]
all_niches_f$spatial_niches <- data_x_y_niches_all$niches
Idents(all_niches_f) <- "spatial_niches"
levels(all_niches_f)
DefaultAssay(all_niches_f) <- 'SCT'
Idents(all_niches_f) <- 'spatial_niches'
DotPlot(all_niches_f, features = rownames(all_niches_f))+RotatedAxis()


pt2 <- table(all_niches_f$named_niches_final, all_niches_f$spatial_niches)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)
pt2$Var2 <- as.double(pt2$Var2)

pt2 <- pt2 %>% filter(Var2 < 6)

cols <- as.data.frame(ArchR::paletteDiscrete(all_niches_f@meta.data[, "named_niches"]))
colnames(cols)<-"colors"

ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)
