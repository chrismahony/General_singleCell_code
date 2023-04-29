library(gsfisher)

getExpressedGenesFromSeuratObject <- function(seurat_object,
                                              clusters,
                                              min.pct=0.1)
{
  expressed <- c()
  for(cluster in clusters)
  {
    # get genes detected in the cluster
    cluster_cells <- names(seurat_object@active.ident[seurat_object@active.ident==cluster])
    clust_pcts <- apply(seurat_object@assays$RNA@data[,cluster_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_clust <- names(clust_pcts[clust_pcts>min.pct])
    
    # get genes detected in the other cells
    other_cells <- names(seurat_object@active.ident[seurat_object@active.ident!=cluster])
    other_pcts <- apply(seurat_object@assays$RNA@data[,other_cells],
                        1, function(x) sum(x>0)/length(x))
    
    detected_in_other_cells <- names(other_pcts[other_pcts>min.pct])
    
    expressed <- c(expressed, detected_in_clust, detected_in_other_cells)
  }
  expressed <- unique(expressed)
  expressed
}

Idents(fibro_only)<-'integrated_snn_res.0.02'
expressed_genes <- getExpressedGenesFromSeuratObject(fibro_only,levels(fibro_only@active.ident), min.pct=0.1)

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
