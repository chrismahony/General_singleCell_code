library(gsfisher)
annotation_gs<-fetchAnnotation(species = "hs")


library(ggfortify)
          library(DESeq2)
          library(dplyr)
          library(ComplexHeatmap)
          library(circlize)
          library(data.table)
          library(gsfisher)
          library(cowplot)

Idents(adult_JIA) <- 'global1'
levels_global <- levels(adult_JIA)

pca_plots <- list()
ht <- list()
go_plot_list()


for (i in 1:length(levels)){

adult_JIA_select<-subset(adult_JIA, idents=levels_global[[i]])

adult_JIA_select$sample_condition<-paste(adult_JIA_select$orig.ident, adult_JIA_select$condition2, sep=".")

cts_fibs<-AggregateExpression(SL_LL_clean, group.by = c("adult_JIA_select"), assays = "RNA", slot = "counts", return.seurat = F)

cts_fibs<-cts_fibs$RNA
cts_fibs<-as.data.frame(cts_fibs)
meta_data=colnames(cts_fibs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all","sample", "condition2")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$disease<-as.factor(meta_data$condition2)

dds <- DESeqDataSetFromMatrix(countData = cts_fibs,
                                  colData = meta_data,
                                  design = ~1)
    
dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

#mingenecount <- quantile(rowSums(counts(dds)), 0.5)
mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

levles <- c("JIA", "Adult")

dds@colData[['condition2']] <- factor(dds@colData[['condition2']],
                                     levels = levels)


design(dds) <- formula(~ condition2)
print(design(dds))
dds <- DESeq(dds, test = "Wald")


#PCA analysis
rld_v3 <- rlog(dds, blind=TRUE)
plotPCA(rld_v3, intgroup="disease")

pcaData <- plotPCA(rld_v3, intgroup=c("sample", "condition2"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plots[[i]] <- ggplot(pcaData, aes(PC1, PC2, color=condition2)) +
  geom_point(aes(shape=condition2, size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

print(pca_plots[[i]])

targetvar <- "condition2"

comps1 <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))
      head(comps1)
      
      ress <- apply(comps1, 1, function(cp) {
        print(cp)
        res <- data.frame(results(dds, contrast=c(targetvar, cp[1], cp[2])))
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

topedges <- 0.05
    ggnet <- Rfast::cora(t(scale_sub_vsd))
    ggnet_gather <- reshape2::melt(ggnet, id.vars = "V1")
    print(dim(ggnet_gather))
    
    ggnet_gather %>%
      filter(value > 0) %>%
      filter(Var1 != Var2) -> ggnet_gather
    
    print(dim(ggnet_gather))
    
    if(nrow(ggnet_gather) > 0) {
      
      edges <- arrange(ggnet_gather, -value)
      print(head(edges))
      print(tail(edges))
      edges <- edges[seq(1, nrow(edges), by = 2), ]
      print(head(edges))
      print(tail(edges))
      qtop <- quantile(edges$value, 1-topedges)
      print(qtop)
      
      edges %>%
        filter(value > qtop) -> top_edges
      
      print(dim(top_edges))
      
      top_edges %>%
        mutate('idx' = paste(Var1, Var2, sep = "_")) -> top_edges# %>%
      #dplyr::distinct(idx, .keep_all = TRUE) -> top_edges
      
      colnames(top_edges) <- c("Source", "Target", "Weight", "Id")
      #top_edges %>%
      #  select(c(Id, Source, Target, Weight)) -> top_edges
      
      # - community detection
      # - igraph definition
      g <- igraph::graph_from_data_frame(top_edges[, c("Source", "Target")],
                                         directed = FALSE)
      g <- igraph::set_edge_attr(g, "weight", value = top_edges$Weight)
      g <- igraph::set_edge_attr(g, "name", value = top_edges$Id)  
      # leiden
      leiden_mod <- igraph::cluster_leiden(g, objective_function = "modularity")
      mods <- data.frame(cbind(igraph::V(g)$name, leiden_mod$membership))
      colnames(mods) <- c("Id", "leiden_gene_cluster")
      imods <- names(table(mods$leiden_gene_cluster)[table(mods$leiden_gene_cluster) > 10])
      print(imods)
      
      if(length(imods) > 0) {
        mods <- filter(mods, leiden_gene_cluster %in% imods)
        head(mods)
        
        row_cl <- data.frame('gene' = mods$Id,
                             'gene_cluster' = paste0("K", mods$leiden_gene_cluster))
        
        row_cl %>%
          arrange(gene_cluster) %>%
          data.frame -> row_cl
        
        rownames(row_cl) <- row_cl$gene
        
              }
    }
    
    rm(g)


colours <- list('sample' = ArchR::paletteDiscrete(adult_JIA@meta.data[, "orig.ident"]),'condition2'= ArchR::paletteDiscrete(adult_JIA@meta.data[, "condition2"]))
col_ann <- HeatmapAnnotation(df = meta_data[, c("sample", "condition2")], col=colours)    
              
    
 if(exists("row_cl")) {
    row_ann <- rowAnnotation(df = row_cl[, -1])    
    
    
   ht[[i]] <-    Heatmap(scale_sub_vsd[rownames(row_cl),], 
              top_annotation = col_ann, 
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = FALSE,
              show_row_names = F,
              show_column_names = F)
   
} 
 
print(ht[[i]])


index <- match(res1$gene, row_cl$gene)
res1$cluster <- row_cl$gene_cluster[index]
res1_cleaned<-na.omit(res1)

index <- match(res1_cleaned$gene, annotation_gs$gene_name)
res1_cleaned$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- unique(res1$gene)
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]
ensemblUni <- na.omit(ensemblUni)
res1_cleaned<-na.omit(res1_cleaned)

go.results <- runGO.all(results=res1_cleaned,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="padj", p_threshold=0.05,
                  species = "hs")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
go_plot[[i]] <- sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

print(go_plot[[i]])

}
