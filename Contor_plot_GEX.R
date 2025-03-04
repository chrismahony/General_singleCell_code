# Extract UMAP coordinates
umap_coords <- seurat_obj@reductions$umap@cell.embeddings

# Extract the gene expression data for the gene of interest
gene_expr <- FetchData(seurat_obj, vars = "Pi16")

# Combine UMAP coordinates and gene expression into a single data frame
umap_data <- data.frame(umap_coords, gene_expr)


# Calculate a filtering threshol (this will need to be customized)
threshold <- mean(umap_data$Pi16)


umap_data$color_expr <- ifelse(umap_data$Pi16 < threshold, NA, umap_data$Pi16)
umap_data_filtered <- umap_data[umap_data$Pi16 >= threshold, ]

# Basic plot with virdis color scale
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = color_expr)) +
  geom_point(alpha = 0.5, size = 1) +  # You can adjust the size and transparency of points if needed
  stat_density_2d(data = umap_data_filtered, aes(fill = ..level..), geom = "polygon", color = "white") +
  scale_fill_viridis_c() +  # Using a viridis color scale for better visualization
  scale_color_viridis_c() +  # Ensure the color scale is the same as the contour fill
  theme_minimal() +
  labs(title = paste("Contour Plot of", "Pi16", "Expression on UMAP"))


# Create the contour plot with a three-color gradient, no gridlines, and a box around the plot
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = color_expr)) +
  geom_point(alpha = 0.5, size = 1) +  # Adjust the size and transparency of points
  # Add density contours with a three-color gradient
  stat_density_2d(data = umap_data_filtered, aes(fill = ..level..), geom = "polygon", color = "white") +
  # Use a three-color gradient for the contour fill
  scale_fill_gradientn(colors = c("blue", "white", "red")) +  # Three-color gradient for contours
  # Use a three-color gradient for the points (from blue to white to red)
  scale_color_gradientn(colors = c("blue", "white", "red")) +  # Three-color gradient for points
  theme_minimal() +
  labs(title = paste("Contour Plot of", "Pi16", "Expression on UMAP")) +
  theme(
    panel.grid.major = element_blank(),   # Remove major gridlines
    panel.grid.minor = element_blank(),   # Remove minor gridlines
    panel.border = element_rect(color = "black", size = 1, fill=NA,),  # Add a black border around the plot
    
  )
