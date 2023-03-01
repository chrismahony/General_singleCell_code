

setwd("your_path") #must contain dirs (one for each sample) with barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz


data.10x = list()
dirs <- list.dirs(".", recursive = FALSE)

for (i in 1:length(dirs)) {
    data.10x[[i]] <- Read10X(data.dir = dirs[[i]])
}
names(data.10x) <- sub('./', '', dirs)
  

scrna.list = list()
samples<-sub('./', '', dirs)

for (i in 1:length(data.10x)) {
    scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=200, project=samples[i]);
    scrna.list[[i]] =NormalizeData(object = scrna.list[[i]]);
    scrna.list[[i]] =ScaleData(object = scrna.list[[i]]);
    scrna.list[[i]] =FindVariableFeatures(object = scrna.list[[i]]);
    scrna.list[[i]] =RunPCA(object = scrna.list[[i]], verbose = FALSE)
    scrna.list[[i]][["percent.mt"]] = PercentageFeatureSet(object=scrna.list[[i]], pattern = "^MT-")
    }
names(scrna.list) <- sub('./', '', dirs)

anchors <- FindIntegrationAnchors(object.list = scrna.list, reduction = "rpca",   dims = 1:50)
aggr <- IntegrateData(anchorset = anchors, dims = 1:50)
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
aggr <- subset(aggr, subset = nFeature_RNA > 500 & percent.mt < 10)
aggr <- FindVariableFeatures(aggr)
aggr <- ScaleData(aggr, verbose = FALSE)
aggr <- RunPCA(aggr, verbose = FALSE)
aggr <- RunUMAP(aggr, dims = 1:50)

aggr <- FindNeighbors(aggr, dims = 1:20)
aggr <- FindClusters(aggr, resolution = c(0.01, 0.05, 0.1, 0.2, 0.3), graph.name = 'integrated_snn')
Idents(aggr)<-'integrated_snn_res.0.01'  
res0.01markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.05'
res0.05markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.1'
res0.1markers<-FindAllMarkers(aggr, only.pos = T)


Idents(aggr)<-'integrated_snn_res.0.2'
res0.2markers<-FindAllMarkers(aggr, only.pos = T)


dev.off()

save.image("you_path/analysis.RData")
