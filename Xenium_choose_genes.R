# Chhose genes that are unique
unique_markers <- straom_markers %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup()


# Pull the raw count matrix from the Seurat object
counts_matrix <- GetAssayData(stroma_clean_h, slot = "counts")

# Subset to include only the unique genes
counts_matrix <- counts_matrix[unique_markers$gene, , drop = FALSE]

# Calculate total counts for each gene and normalize by the total number of cells
total_counts <- rowSums(counts_subset)  
num_cells <- ncol(counts_subset)       
mean_expression <- total_counts / num_cells

# Create a DataFrame for easy viewing
mean_expression_df <- data.frame(
  gene = rownames(counts_subset),
  mean_expression = mean_expression
)

df <- mean_expression_df %>% filter(mean_expression > 0.1 & mean_expression < 100)
df
