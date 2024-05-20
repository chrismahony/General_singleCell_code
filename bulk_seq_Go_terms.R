
counts2 <- counts_all %>% select(-c(Chr, Start,End, Strand, Length))


counts2_F <- counts2 %>% select(c(1:2, 4:5, 7:8, 10:12))
colnames(counts2_F) <- c("R1A1", "R1A2", "R1C1", "R1C2", "EV1", "EV2", "EV4", "R1A4", "R1C4")

counts2_F <- counts2_F[c(5,6,7,1,2,8,3,4,9)]

meta_data=colnames(counts2_F)
meta_data<-as.data.frame(meta_data)

condition <- c(rep("EV", 3), rep("R1A", 3), rep("R1C", 3))
meta_data$condition <- condition
colnames(meta_data)<-c("sample", "condition")

meta_data$sample<-as.factor(meta_data$sample)
meta_data$disease<-as.factor(meta_data$condition)

library(sva)
batch <- c(1,1,2,1,1,2,1,1,2)


adjusted <- ComBat_seq(counts2_F, batch=batch, group=NULL)

dds <- DESeqDataSetFromMatrix(countData = adjusted,
                                  colData = meta_data,
                                  design = ~1)

print(quantile(rowSums(counts(dds))))

mingenecount <- quantile(rowSums(counts(dds)), 0.5)
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

dds@colData[['condition']] <- factor(dds@colData[['condition']])


design(dds) <- formula(~ condition)
print(design(dds))
dds <- DESeq(dds, test = "Wald")
plotDispEsts(dds)

#PCA analysis
rld_v3 <- rlog(dds, blind=TRUE)
plotPCA(rld_v3, intgroup="condition")

pcaData <- plotPCA(rld_v3, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(aes(shape=condition, size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


rld_mat_v3 <- assay(rld_v3)    
rld_cor_v3 <- cor(rld_mat_v3)
library(pheatmap)
pheatmap(rld_cor_v3)

targetvar <- "condition"

comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(DESeq2::results(dds, contrast=c(targetvar, cp[1], cp[2])))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- paste0(cp[1], "_vs_", cp[2])
        res
      })
      
      res1 <- Reduce(rbind, ress)

res1 %>% 
      filter(padj < 0.05) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

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

col_ann <- HeatmapAnnotation(df = meta_data)
Heatmap(scale_sub_vsd, 
              top_annotation = col_ann, 
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = F)





### filter DEGs as you want###

res %>% 
  filter(padj < 0.01) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> subres


# add 'cluster' column to DEGs ###

subres$cluster <- "NO"
subres$cluster[DEGs$log2FoldChange > 0.001] <- "UP"
subres$cluster[DEGs$log2FoldChange < -0.001] <- "DOWN"


library(gsfisher)
#annotation_gs<-fetchAnnotation(species = "hs")

index <- match(subres$gene, annotation_gs$gene_name)
subres$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- unique(subres$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)

go.results <- runGO.all(results=subres,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)


