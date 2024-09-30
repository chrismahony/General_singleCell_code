library(harmony)

seuratObj <- RunHarmony(seuratObj, group.by.vars=c("meta.data.slot1", "meta.data.slot2", "meta.data.slot3", ...))
seuratObj <- RunUMAP(seuratObj, reduction = "harmony")
DimPlot(seuratObj)
seuratObj <- FindNeighbours(seuratObj, dims=1:30, reduction = "harmony")
seuratObj <- FindClusters(seuratObj, res=c(0.05, 0.1, 0.2,0.3,...)
