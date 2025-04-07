library(fgsea)
library(data.table)
library(ggplot2)
data(examplePathways)
library(org.Hs.eg.db)
data(exampleRanks)

gene_ranks1 <- deg_df1$stat
names(gene_ranks1) <- deg_df1$gene
gene_ranks1 <- sort(gene_ranks1, decreasing = TRUE)
fgsea1 <- fgsea(pathways=examplePathways, stats = gene_ranks1)

# Simulated data, but this should be DEG lists. I.e. Upregulayted genes in treated1 vs control etc.
ranks1 <- sample(exampleRanks, size=1000)
ranks2 <- sample(exampleRanks, size=1000)
ranks3 <- sample(exampleRanks, size=1000)


fgsea1 <- fgsea(pathways=examplePathways, stats = ranks1)
fgsea2 <- fgsea(pathways=examplePathways, stats = ranks2)
fgsea3 <- fgsea(pathways=examplePathways, stats = ranks3)


get_nes_vector <- function(fgsea_result, sample_name) {
  fgsea_result %>%
    dplyr::select(pathway, NES) %>%
    dplyr::rename(!!sample_name := NES)
}

nes_list <- list(
  get_nes_vector(fgsea1, "Sample1"),
  get_nes_vector(fgsea2, "Sample2"),
  get_nes_vector(fgsea3, "Sample3")
)

library(tidyverse)

nes_mat <- reduce(nes_list, full_join, by = "pathway") %>%
  column_to_rownames("pathway")  

nes_mat

nes_mat[is.na(nes_mat)] <- 0

col_fun <- colorRamp2(
  c(min(nes_mat), 0, max(nes_mat[nes_mat > 0])),  
  c("white", "white", "red")
)

Heatmap(nes_mat, show_row_names = F, border=T)

Heatmap(
  nes_mat,
  name = "NES",
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = T,
  border=T
)
