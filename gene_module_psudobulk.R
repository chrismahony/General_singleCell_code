library(tidyverse)
#create DESeq2 object first
myel_hs$cluster_type<-paste(myel_hs$named, myel_hs$TYPE, sep=".")
cts_all<-AggregateExpression(myel_hs, group.by = c("cluster_type"), assays = "RNA", slot = "counts", return.seurat = F)

#get meta_data read for dseq2
cts_all<-cts_all$RNA
#cts.t<-t(cts)
cts_all<-as.data.frame(cts_all)
meta_data=colnames(cts_all)
meta_data<-as.data.frame(meta_data)
library(splitstackshape)
meta_data$to_split<-meta_data$meta_data
meta_data<-cSplit(meta_data, splitCols = "to_split", sep=".")
colnames(meta_data)<-c("all_details", "cluster", "condiiton")
#meta_data$sample_global<-paste(meta_data$sample, meta_data$global, sep="_")

library(DESeq2)
dds_all<-DESeqDataSetFromMatrix(countData = cts_all, colData=meta_data, design = ~ condiiton)
keep<-rowSums(counts(dds_all))>=100
dds_all<-dds_all[keep,]
dds_all <- DESeq(dds_all)
normalized_counts <- counts(dds_all, normalized=TRUE)

resultsNames(dds_all)
deseq2Results <- results(dds_all, contrast=c("SF", "Tissue", "Blood"))


library(tidyverse)
normalized_counts_n<-as.data.frame(normalized_counts)
normalized_counts_n$gene_ID<-rownames(normalized_counts_n)
Exp_table_long <- normalized_counts_n %>% 
  pivot_longer(cols = !gene_ID, names_to = "library", values_to = "nomalized_counts") %>% 
  mutate(logCounts = log10(nomalized_counts + 1))

Exp_table_long_averaged <- Exp_table_long %>% 
   group_by(gene_ID, library) %>% 
  summarise(mean.logCounts = mean(logCounts)) %>% 
  ungroup()  

Exp_table_long_z <- Exp_table_long_averaged %>% 
  group_by(gene_ID) %>% 
  mutate(z.score = (mean.logCounts - mean(mean.logCounts))/sd(mean.logCounts)) %>% 
  ungroup()

high_var_genes_no_average <- Exp_table_long_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logCounts)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes_no_average)
dim(high_var_genes_no_average)

#choose number of highly varible genes
high_var_genes10000 <- high_var_genes_no_average %>% 
  slice_max(order_by = var, n = 10000)

Exp_table_long_z_high_var <- Exp_table_long_z %>% 
  filter(gene_ID %in% high_var_genes10000$gene_ID)

z_score_wide <- Exp_table_long_z_high_var %>% 
  select(gene_ID, library, z.score) %>% 
  pivot_wider(names_from = library, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
head(z_score_wide)

cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)


number_of_tissue_stage <- ncol(z_score_wide) - 1
number_of_tissue_stage

cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA

edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_matrix)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_of_tissue_stage-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_tissue_stage-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 

head(edge_table)

edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )  
  
  
 #filter edges and choose cluster resolution
library(tidyverse)
library(igraph)

edge_table_select <- edge_table %>% 
  filter(r >= 0.4)


my_network <- graph_from_data_frame(
  edge_table_select,
    directed = F
)


modules <- cluster_leiden(my_network, resolution_parameter = 1.25, 
                          objective_function = "modularity")


#check length and re run clus res if needed
length(modules)


my_network_modules <- data.frame(
  gene_ID = names(membership(modules)),
  module = as.vector(membership(modules)) 
) 



Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_z_high_var %>% 
  inner_join(my_network_modules, by = "gene_ID")

modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, library) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

library(splitstackshape)
modules_mean_z$sample<-modules_mean_z$library
modules_mean_z<-cSplit(modules_mean_z, splitCols="library", sep=".")

module_peak_exp <- modules_mean_z %>% 
  group_by(module, library_2) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp


modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, 0.95)

modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 3 ~ 3,
    mean.z < -3 ~ -3,
    T ~ mean.z
  ))


modules_mean_z_wide <- modules_mean_z %>% 
  select(module, sample, mean.z.clipped
) %>% 
  pivot_wider(names_from = sample, values_from = mean.z.clipped
)

modules_mean_z_wide$module<-NULL



library(ComplexHeatmap)

library(circlize)
#add heatmap annotaitons (if needed)
ha = HeatmapAnnotation(
    type = c(rep("Tissue", 9), rep("SF", 9), rep("Blood", 9)  ), 
    sample = 1:27,
    col = list(type = c("Tissue" = "red", "SF" = "green", "Blood" = "blue"))
)


Heatmap(modules_mean_z_wide, cluster_columns = F,  top_annotation=ha)
  
  
