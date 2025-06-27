#PBMC1 is you main object with all cells
#objs will be you rseurat objs

objs <- list()
for (i in 1:length(umap_embeddings)){

common_cells <- intersect(colnames(PBMC1), rownames(umap_embeddings[[i]]))

objs[[i]] <- subset(PBMC1, cells = common_cells)

umap_embedding_f <- umap_embeddings[[i]][common_cells, ]

objs[[i]][["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap_embedding_f), 
                                              key = "UMAP_", assay = DefaultAssay(objs[[i]]))

objs[[i]] <- objs[[i]] %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()

DimPlot(objs[[i]], raster=FALSE)

}

names(objs) <- names(umap_embeddings)

objs[["stromal"]] %>% DimPlot()

objs[["stromal"]] %>% DimPlot(group.by="subclusters")
