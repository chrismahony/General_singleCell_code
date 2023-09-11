library(ggfortify)
          library(DESeq2)
          library(dplyr)
          library(ComplexHeatmap)
          library(circlize)
          library(data.table)
          library(gsfisher)
          library(cowplot)


#select cells you want
SL_LL_clean<-subset(SL_LL, idents=c("Cthrc1_C1qtnf3", "Cxcl5_Mmp3" ,    "Prg4_Clic5" ,    "Gsn_Plpp3"    ,"Cxcl10_Ccl7"  ,  "Cd34_Pi16"   ,   "Smoc2_Ccl11" ,  "Eln_Gdf10"  ,                      "Sfrp5_Cldn1"   ))


#pseudobulk and create meta.data
SL_LL_clean$sample_inflammation_model<-paste(SL_LL_clean$orig.ident, SL_LL_clean$inflammation, SL_LL_clean$model, sep=".")
Idents(SL_LL_clean)<-'sample_inflammation_model'
levels(SL_LL_clean)

cts_fibs<-AggregateExpression(SL_LL_clean, group.by = c("sample_inflammation_model"), assays = "RNA", slot = "counts", return.seurat = F)

cts_fibs<-cts_fibs$RNA
cts_fibs<-as.data.frame(cts_fibs)
meta_data=colnames(cts_fibs)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all","sample", "disease", "model")
meta_data$all<-as.factor(meta_data$all)
meta_data$sample<-as.factor(meta_data$sample)
meta_data$disease<-as.factor(meta_data$disease)
meta_data$model<-as.factor(meta_data$model)

#meta_CIA<-meta_data[meta_data$model == "CIA",]
#cts_fibs_CIA<-cts_fibs[,meta_CIA$all]

#meta_AIA<-meta_data[meta_data$model == "AIA",]
#cts_fibs_AIA<-cts_fibs[,meta_AIA$all]

#meta_STIA<-meta_data[meta_data$model == "STIA",]
#cts_fibs_STIA<-cts_fibs[,meta_STIA$all]

meta_peak<-meta_data[meta_data$disease == "peak",]
meta_rest<-meta_data[meta_data$disease == "rest",]
meta_resing<-meta_data[meta_data$disease == "resing",]

meta_peak<-rbind(meta_peak,meta_rest)
meta_resing<-rbind(meta_resing,meta_rest)

cts_fibs_peak<-cts_fibs[,meta_peak$all]
cts_fibs_resing<-cts_fibs[,meta_resing$all]

#start of analysis
meta<-meta_resing
cts<-cts_fibs_resing
#levels<-c("rest","initiation","peak","resing","resed","persis")
#levels<-c("rest","peak")
levels<-c("rest","resing")
#levels<-c("rest","peak","Eresing","resing","resed")

dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = meta,
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

dds@colData[['disease']] <- factor(dds@colData[['disease']],
                                     levels = levels)

dds@colData[['disease']] <- as.factor(dds@colData[['disease']])

design(dds) <- formula(~ disease)
print(design(dds))
dds <- DESeq(dds, test = "Wald")


#PCA analysis
rld_v3 <- rlog(dds, blind=TRUE)
plotPCA(rld_v3, intgroup="disease")

pcaData <- plotPCA(rld_v3, intgroup=c("disease", "model"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=disease)) +
  geom_point(aes(shape=model, size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(pcaData, aes(PC1, PC2, color=model)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


rld_mat_v3 <- assay(rld_v3)    
rld_cor_v3 <- cor(rld_mat_v3)
library(pheatmap)
pheatmap(rld_cor_v3)


    
targetvar <- "disease"

comps1 <- data.frame(t(combn(unique(as.character(meta[[targetvar]])), 2)))
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
      filter(padj < 0.01) %>%
      mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
      arrange(desc(abs(score))) -> subres

#optional
#subres_up<-subres[subres$log2FoldChange > 1,]
#subres_down<-subres[subres$log2FoldChange < -1,]
#subres<-rbind(subres_up, subres_down)

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

ss_sm <- meta_data[, c("model", "disease", "sample")]
col_ann <- HeatmapAnnotation(df = ss_sm)

topedges <- 0.05
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
#AIA    
#   col.order<-c("rest_neg_renamed.rest.AIA", "Day0neg_cr_V5.rest.AIA", #"AIA_rest1_cd45neg.rest.AIA","AIA_infla1_cd45neg.peak.AIA", #"AIA_infla2_cd45neg.peak.AIA","Day2_neg3.peak.AIA" , "Day4_neg1.Eresing.AIA" , #"Day4_neg3.Eresing.AIA","Day4_neg2.Eresing.AIA",  "Day7_neg1.resing.AIA"  #,"Day7_neg2.resing.AIA", "Day7_neg3.resing.AIA" , "D14_1_neg_renamed.resed.AIA", #"D14_2_neg_renamed.resed.AIA", "D14_3_neg_renamed.resed.AIA"  )
#    scale_sub_vst_mat<-scale_sub_vst_mat[ , col.order]
    
#    colnames<-as.data.frame(colnames(scale_sub_vst_mat))
#    colnames<-cSplit(colnames, splitCols = "colnames(scale_sub_vst_mat)", sep=".")
#  colnames(colnames)<-c("sample", "disease", "model")
    
#STIA   
 col.order<-c("CD45N.rest.STIA", "S8_CONTROLnegB_r.rest.STIA", "S8_ControlnegC.rest.STIA", "S2_DAY1negA_cmAB2.initiation.STIA", "S4_DAY1negB_cmAB2.initiation.STIA", "S6_DAY1negC_cmAB2.initiation.STIA", "S2_DAY8negA.peak.STIA",  "S4_Day8negB.peak.STIA", "S6_DAY8negC.peak.STIA", "S2_Day15negB.resing.STIA","S4_Day15negC.resing.STIA","S8_DAY15negA.resing.STIA" , "S2_Day22negA.resed.STIA",  "S4_Day22negB.resed.STIA" , "S6_Day22negC.resed.STIA", "S2_Day28negC.persis.STIA", "S6_Day28negA.persis.STIA", "S8_Day28negB.persis.STIA"
     
   )
 
 
  #peak comparison   
 col.order<-c(  "S8_CONTROLnegB_r.rest.STIA" , "S8_ControlnegC.rest.STIA",  "CD45N.rest.STIA", "Con1_CD45neg.rest.CIA",      
"Con2_CD45neg.rest.CIA" ,      "Con3_CD45neg.rest.CIA",       "Con4_CD45neg.rest.CIA" ,     "Day0neg_cr_V5.rest.AIA"   , "AIA_rest1_cd45neg.rest.AIA",  "rest_neg_renamed.rest.AIA"   , "S2_DAY8negA.peak.STIA" ,      "S4_Day8negB.peak.STIA"     , 
 "S6_DAY8negC.peak.STIA",        "Infla1_CD45neg.peak.CIA"  , "Infla2_CD45neg.peak.CIA"  ,   "Infla3_CD45neg.peak.CIA" , "AIA_infla1_cd45neg.peak.AIA" ,"AIA_infla2_cd45neg.peak.AIA", "Day2_neg3.peak.AIA"   )
 
 #scale_sub_vst_mat<-scale_sub_vst_mat[ , col.order]

          

    
    colnames<-as.data.frame(colnames(scale_sub_vst_mat))
    colnames<-cSplit(colnames, splitCols = "colnames(scale_sub_vst_mat)", sep=".")
  colnames(colnames)<-c("sample", "disease", "model")    

 ss_sm <- colnames[, c("model", "disease", "sample")]
col_ann <- HeatmapAnnotation(df = ss_sm)


 if(exists("row_cl")) {
    row_ann <- rowAnnotation(df = row_cl[, -1])    
    
    draw(
      Heatmap(scale_sub_vst_mat[rownames(row_cl),], 
              top_annotation = col_ann, 
              right_annotation = row_ann,
              col=colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 4), 
              cluster_columns = F,
              cluster_rows = FALSE,
              show_row_names = F,
              show_column_names = F))
} 
 

#output_DEG_CIAs=res1
#output_mods_CIA=mods

#output_DEG_AIAs=res1
#output_mods_AIA=mods

#output_DEG_STIAs=res1
#output_mods_STIA=mods

#write.csv(output_DEG_CIAs, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_DEG_CIAs.csv")
#write.csv(output_mods_CIA, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_mods_CIAs.csv")

#write.csv(output_DEG_STIAs, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_DEG_STIAs.csv")
#write.csv(output_mods_STIA, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_mods_STIAs.csv")

#write.csv(output_DEG_AIAs, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_DEG_AIAs.csv")
#write.csv(output_mods_AIA, "/rds/projects/c/croftap-sitia-cite-seq-tc/atlas_pipeline_v2/output_mods_AIAs.csv")


library(gsfisher)
annotation_gs<-fetchAnnotation(species = "mm")


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
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=5, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)

