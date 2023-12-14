
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

