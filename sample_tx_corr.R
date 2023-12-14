
########pb as normal##########

stia2021_rna$sample_condition_cluster<-paste(stia2021_rna$sample_id, stia2021_rna$condition,stia2021_rna$cluster.name, sep=".")
stia2021_rna$condition_cluster<-paste(stia2021_rna$condition,stia2021_rna$cluster.name, sep=".")

cts_new<-AggregateExpression(stia2021_rna, group.by = c("condition_cluster"), assays = "RNA", slot = "counts", return.seurat = F)

rm(list=ls()[! ls() %in% c("cts_new","stia2021_rna")])
gc()

cts_new<-cts_new$RNA
cts_new<-as.data.frame(cts_new)
meta_data_new=colnames(cts_new)
meta_data_new<-as.data.frame(meta_data_new)
library(splitstackshape)
meta_data_new$to_split<-meta_data_new$meta_data_new
meta_data_new<-cSplit(meta_data_new, splitCols = "to_split", sep=".")
colnames(meta_data_new)<-c("all", "condition", "cluster")
meta_data_new$all<-as.factor(meta_data_new$all)
meta_data_new$sample<-as.factor(meta_data_new$sample)
meta_data_new$cluster<-as.factor(meta_data_new$cluster)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts_new,
                                  colData = meta_data_new,
                                  design = ~1)

dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

mingenecount <- 200
#maxgenecount <- quantile(rowSums(counts(dds)), 0.99)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount #& rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['condition']] <- as.factor(dds@colData[['condition']])

design(dds) <- as.formula(paste0("~", "condition"))

print(design(dds))

dds <- DESeq(dds, test = "Wald")

meta <- meta_data_new
print(resultsNames(dds))
targetvar <- "condition"
comps <- data.frame(t(combn(unique(as.character(meta[[targetvar]])), 2)))
head(comps)

ress <- apply(comps, 1, function(cp) {
  print(cp)
  res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
  res[["gene"]] <- rownames(res)
  res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
  res
})



res <- Reduce(rbind, ress)

head(res)
comps

res %>% 
  filter(padj < 0.01) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres

head(subres)
dim(subres)

length(unique(subres$gene))


library(ComplexHeatmap)

      if(length(unique(subres$gene)) > 10) {
      vsd <- tryCatch({
        vst(dds, blind=TRUE)
      }, error=function(e) {
        message(e)
        print(e)
        return(NULL)
      })
      
      if(!is.null(vsd)) {
        print(dim(assay(vsd)))
        print(head(assay(vsd), 3))
        vsd_mat <- assay(vsd)
        
        feats <- unique(subres$gene)
        print(length(feats))
        
        # Sub-set matrix to relevant features
        sub_vsd_mat <- vsd_mat[rownames(vsd_mat) %in% feats, ]
        scale_sub_vsd <- t(scale(t(sub_vsd_mat)))
        head(scale_sub_vsd)
        dim(scale_sub_vsd)
      }
      }


######### PC analysis for each sample based on GEX ########

pc <- prcomp(t(scale_sub_vsd),
             center = TRUE,
            scale. = TRUE)


m <- pc$x %>% as.matrix() 


library(harmony)
harmony_embeddings <- harmony::HarmonyMatrix(
    m, meta_data_new, c('condition'), do_pca = F, verbose= TRUE
)

#####Heatmap   ######

colours <- list('cluster' = ArchR::paletteDiscrete(stia2021_rna@meta.data[, "cluster.name"]),'condition'= ArchR::paletteDiscrete(stia2021_rna@meta.data[, "condition"]))
col_ann <- HeatmapAnnotation(df = meta_data_new[,c("cluster", "condition")], col=colours)

Heatmap(cor(t(harmony_embeddings)),  show_row_names = F,
        show_column_names = F, column_dend_reorder = TRUE, 
        name="harmony space\ncorrelation",
        raster_quality = 3,
        use_raster = TRUE,top_annotation = col_ann )
