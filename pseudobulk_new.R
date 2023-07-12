library(ggfortify)
          library(DESeq2)
          library(dplyr)
          library(ComplexHeatmap)
          library(circlize)
          library(data.table)
          library(gsfisher)
          library(cowplot)


merged_down$orig.ident_named_data<-paste(merged_down$orig.ident, merged_down$named, merged_down$data, sep=".")
table(merged_down$orig.ident_named_data)

cts<-AggregateExpression(merged_down, group.by = c("orig.ident_named_data"), assays = "RNA", slot = "counts", return.seurat = F)

cts<-cts$RNA
cts<-as.data.frame(cts)
meta_data=colnames(cts)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all", "sample", "cluster", "data")
meta_data$cluster_data<-paste(meta_data$cluster, meta_data$data, sep="_")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$cluster<-as.factor(meta_data$cluster)
meta_data$data<-as.factor(meta_data$data)
meta_data$cluster_data<-as.factor(meta_data$cluster_data)


dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = meta_data,
                                  design = ~1)
    
dds <- scran::computeSumFactors(dds)
print(dds)
print(quantile(rowSums(counts(dds))))

mingenecount <- 200
maxgenecount <- quantile(rowSums(counts(dds)), 0.999)
dim(counts(dds))
# Subset low-expressed genes
keep <- rowSums(counts(dds)) > mingenecount & rowSums(counts(dds)) < maxgenecount
dds <- dds[keep, ]
print(quantile(rowSums(counts(dds))))
dim(dds)

design(dds) <- formula(~ cluster_data)
print(design(dds))
dds <- DESeq(dds, test = "Wald")
    
print(resultsNames(dds))

targetvar <- "cluster_data"

comps <- resultsNames(dds)
      ress <- lapply(comps[-1], function(cp) {
        print(cp)
        res <- data.frame(results(dds, name=cp))
        res[["gene"]] <- rownames(res)
        res[["comparison"]] <- cp
        res
      })

      
res <- Reduce(rbind, ress)
    
head(res)      

res %>% 
      filter(padj < 0.01) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

subres_up<-subres[subres$log2FoldChange > 1,]
subres_down<-subres[subres$log2FoldChange < -1,]
subres<-rbind(subres_up, subres_down)

head(subres)
dim(subres)
length(unique(subres$gene))

library(ComplexHeatmap)

#use vsd for <10 uniqe DEGs

deseq2VST <- vst(dds, blind=T)

feats <- unique(subres$gene)
print(length(feats))
        
# Sub-set matrix to relevant features
deseq2VST <- assay(deseq2VST)
deseq2VST<-as.matrix(deseq2VST)
sub_vst_mat <- deseq2VST[rownames(deseq2VST) %in% feats, ]
scale_sub_vst_mat <- t(scale(t(sub_vst_mat)))
head(scale_sub_vst_mat)
dim(scale_sub_vst_mat)


ss_sm <- meta_data[, c("sample", "data", "cluster")]
col_ann <- HeatmapAnnotation(df = ss_sm)


nrow(scale_sub_vst_mat)
        
draw(Heatmap(scale_sub_vst_mat, top_annotation = col_ann,
             col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                  row_names_gp = gpar(fontsize = 4), 
                  cluster_columns = F, use_raster=F))


 for(v in colnames(ss_sm)) {
          plot(
            autoplot(prcomp(t(scale_sub_vst_mat)), 
                     data=ss_sm,
                     colour = v,
                     label = FALSE, 
                     label.size = 3) +
              theme_classic() +
              theme(legend.position = "bottom")
          )}

        

    ggnet <- Rfast::cora(t(scale_sub_vst_mat))
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
      
      print(head(top_edges))
      print(dim(top_edges))
      
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
  
    
if(exists("row_cl")) {
    row_ann <- rowAnnotation(df = row_cl[, -1])    
    
    draw(
      Heatmap(scale_sub_vst_mat[rownames(row_cl),], 
              top_annotation = col_ann, 
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = T,
              cluster_rows = FALSE,
              show_row_names = F,
              show_column_names = F))
}
