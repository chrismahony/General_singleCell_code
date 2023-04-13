setwd("/rds/projects/m/mcgetthm-rnaseq-arthritis/")
dirs <-dir("/rds/projects/m/mcgetthm-rnaseq-arthritis/", pattern = "S")
dirs <- dirs[dirs != c("FASTQ", "FASTQs_CM")]
samples<- dirs
dirs <- paste0( "./", dirs, "/outs/") 
data.10x = list()

# percentage of doublets from https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

db_pct<-as.numeric(c("0.03", "0.08"))


#to configure doublet finder

#create seurat list and QC
for (i in 1:2) {
data.10x[[i]] = load10X(dataDir =dirs[[i]])
data.10x[[i]] = autoEstCont(data.10x[[i]])
data.10x[[i]] = adjustCounts(data.10x[[i]], roundToInt=T )
data.10x[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=3, min.features=0, project=samples[i]);
data.10x[[i]][["percent.mt"]] = PercentageFeatureSet(object=data.10x[[i]], pattern = "^mt-");
data.10x[[i]]<-subset(x = data.10x[[i]], subset = nFeature_RNA > 200 & nFeature_RNA <7000 & percent.mt < 10);
data.10x[[i]] =NormalizeData(object = data.10x[[i]]);
data.10x[[i]] =ScaleData(object = data.10x[[i]]);
data.10x[[i]] =FindVariableFeatures(object = data.10x[[i]]);
data.10x[[i]] =RunPCA(object = data.10x[[i]], verbose = FALSE);
data.10x[[i]] =RunUMAP(object = data.10x[[i]], dims=1:30)
}

#calclate optimal pk
sweep.stats.list <- list()
optimal.pk.list<-list()
for (i in 1:length(data.10x)) {
  seu_temp <- data.10x[[i]]
  sweep.res.list <- paramSweep_v3(seu_temp, PCs = 
  1:30)
  sweep.stats <- summarizeSweep(sweep.res.list)
  bcmvn <- find.pK(sweep.stats)
  sweep.stats.list[[i]] <- sweep.stats
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  optimal.pk.list[[i]] <- optimal.pk
}

#claisfy doublets
optimal.pk.list_s<-optimal.pk.list
#next to run
db_pct<-as.numeric(c("0.03", "0.08"))
optimal.pk.list<-as.numeric(optimal.pk.list)
cleaned.list<-list()
for (i in 1:length(data.10x)) {
  seu_temp <- data.10x[[i]]
  nExp_poi <-  db_pct[[i]]*nrow(seu_temp@meta.data)
  seu_temp <- doubletFinder_v3(seu_temp, PCs = 
  seu_temp@commands$RunUMAP.RNA.pca$dims, pN = 0.25, pK = optimal.pk.list[i], nExp = nExp_poi, reuse.pANN = FALSE)
  cleaned.list[[i]] <- seu_temp
}

names(cleaned.list)<-samples


#remove doublets and re process data
for (i in 1:length(cleaned.list)) {
DF.name = colnames(cleaned.list[[i]]@meta.data)[grepl("DF.classification", colnames(cleaned.list[[i]]@meta.data))];
cleaned.list[[i]] = cleaned.list[[i]][, cleaned.list[[i]]@meta.data[, DF.name] == "Singlet"];
cleaned.list[[i]] =ScaleData(object = cleaned.list[[i]]);
cleaned.list[[i]] =FindVariableFeatures(object = cleaned.list[[i]]);
cleaned.list[[i]] =RunPCA(object = cleaned.list[[i]], verbose = FALSE)
}
