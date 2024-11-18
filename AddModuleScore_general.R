library(Seurat)

#filter df from bulk RNAseq DEGs

df_f <- subres %>% filter(log2FoldChange > 1 & padj < 0.05)

#if multiple comparisons have been made
df_f <- subres %>% filter(comparison == "Dex_vs_LPS" & log2FoldChange > 1 & padj < 0.05)  #adulst cot off as required, you will need to try different ones

### !! make sure you take the correct log2FOldChnage direcetion. log2FoldChange > 1 in Dex_vs_LPS are gene increased in Dex ###

#Add to Seurat obj

obj <- AddModuleScore(obj, features = list(df_f$gene_name), name="module")


#Visulise

DotPlot(obj, features="module1")  #Add 1 to name
