options(bitmapType='cairo')
library(Seurat)


# Generate transfer anchors
transfer.anchors <- FindTransferAnchors(reference = obj.rna, query = stia2021_rna_s, features = VariableFeatures(object = obj.rna),
    reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata =obj.rna$Cluster_name,
    weight.reduction = stia2021_rna_s[["pca"]], dims = 1:30)

# Add this to Seurat obj
stia2021_rna_s <- AddMetaData(stia2021_rna_s, metadata = celltype.predictions)


# Calculate mapping proportions
predictions <- table(stia2021_rna_s$cluster.name, stia2021_rna_s$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p3 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation XXX") + ylab("Predicted cell type label YYY") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p3


# optional- Generate Heatmap
predictions2<-predictions %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 

library(tidyverse)
predictions2 <-predictions2 %>% remove_rownames %>% column_to_rownames(var="Var1")


Heatmap(predictions2, cluster_rows = T, colorRamp2(c(-max(predictions2), 0, max(predictions2)), c('white', 'white', 'black')), cluster_columns= T)
