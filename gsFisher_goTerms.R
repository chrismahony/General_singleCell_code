library(gsfisher)

expressed_genes<-rownames(fibro_only)
annotation_gs <- fetchAnnotation(species="mm", ensembl_version=NULL, ensembl_host=NULL)

index <- match(markers$gene, annotation_gs$gene_name)
markers$ensembl <- annotation_gs$ensembl_id[index]

FilteredGeneID <- expressed_genes
index <- match(FilteredGeneID, annotation_gs$gene_name)
ensemblUni <- annotation_gs$ensembl_id[index]

seurat_obj.res <- markers
seurat_obj <- fibro_only
seurat_obj.res <- seurat_obj.res[!is.na(seurat_obj.res$ensembl),]
ensemblUni <- na.omit(ensemblUni)

go.results <- runGO.all(results=seurat_obj.res,
                  background_ids = ensemblUni, gene_id_col="ensembl", gene_id_type="ensembl", sample_col="cluster", p_col="p_val_adj", p_threshold=0.05,
                  species = "mm")
go.results <- filterGenesets(go.results)
go.results.top <- go.results %>% group_by(cluster) %>% top_n(n=32, -p.val)
sampleEnrichmentDotplot(go.results.top, selection_col = "description", selected_genesets = unique(go.results.top$description), sample_id_col = "cluster", fill_var = "odds.ratio", maxl=50, title="Go term",rotate_sample_labels = T)
