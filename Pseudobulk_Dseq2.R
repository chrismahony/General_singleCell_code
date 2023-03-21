#paste together all the varibles you are interested in
v3$global_sample_condition<-paste(v3$global , v3$sample, v3$condition, sep=".")

#AggregateExpression will sum counts when slot is set to "counts"
cts_v3<-AggregateExpression(v3, group.by = c("global_sample_condition"), assays = "RNA", slot = "counts", return.seurat = F)

#get meta_data read for dseq2
cts_v3<-cts_v3$RNA
#cts.t<-t(cts)
cts_v3<-as.data.frame(cts_v3)
meta_data=colnames(cts_v3)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all_deatils", "global", "sample", "condition")
meta_data$sample_global<-paste(meta_data$sample, meta_data$global, sep="_")


dds_v3<-DESeqDataSetFromMatrix(countData = cts_v3, colData=meta_data, design = ~ condition)
keep<-rowSums(counts(dds_v3))>=10
dds_v3<-dds_v3[keep,]


#PCA analysis
rld_v3 <- rlog(dds_v3, blind=TRUE)
plotPCA(rld_v3, intgroup="sample_global")

pcaData <- plotPCA(rld_v3, intgroup=c("meta_data", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


rld_mat_v3 <- assay(rld_v3)    
rld_cor_v3 <- cor(rld_mat_v3)
library(pheatmap)
pheatmap(rld_cor_v3)

#dseq2
dds_v3 <- DESeq(dds_v3)
deseq2Results <- results(dds_v3)
deseq2ResDF <- as.data.frame(deseq2Results)
summary(deseq2Results)
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)
#strict
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .01 & abs(deseq2ResDF$log2FoldChange) > 3,])


deseq2VST <- vst(dds_v3)
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

library(ComplexHeatmap)
Heatmap(deseq2VST,cluster_columns =F)
