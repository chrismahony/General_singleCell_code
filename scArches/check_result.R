#read in after scArches
meta_amp2_scarhces_labels_df <- read_csv("/rds/projects/c/croftap-celldive01/amp2/prepare_scArches/meta_amp2_scarhces_labels_df.csv")


meta_amp2_scarhces_labels_df$cell_type_pred %>% table()

meta_amp2_scarhces_labels_df %>% colnames()


predictions <- table(meta_amp2_scarhces_labels_df$cell_type_pred, meta_amp2_scarhces_labels_df$cluster_name)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


df <- predictions %>% pivot_wider(names_from = Var1, values_from = Freq) %>% as.data.frame()
rownames(df) <- df$Var2
df$Var2 <- NULL
df %>% Heatmap()



library("scProportionTest")

Idents(fibs_amp2)<-'sample'
CTAP_donor_mapping <- read_excel("/rds/projects/c/croftap-celldive01/amp2/processed_output_04-11-2023/CTAP_donor_mapping.xlsx")

donors <- CTAP_donor_mapping$donor
donors[-43]

fibs_amp2_RA<-subset(fibs_amp2, idents=donors[-43])

table(fibs_amp2$sample)

fibs_amp2_RA <- fibs_amp2_RA %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

meta_amp2_scarhces_labels_df$...1 <- gsub('.{2}$', '', meta_amp2_scarhces_labels_df$...1)
meta_amp2_scarhces_labels_df$...1 %>% tail()

meta_amp2_scarhces_labels_df_f <- meta_amp2_scarhces_labels_df[meta_amp2_scarhces_labels_df$cell %in% colnames(fibs_amp2_RA),]

nrow(meta_amp2_scarhces_labels_df_f)
ncol(fibs_amp2_RA)



predictions <- table(meta_amp2_scarhces_labels_df_f$cell_type_pred, meta_amp2_scarhces_labels_df_f$cluster_name)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


df <- predictions %>% pivot_wider(names_from = Var1, values_from = Freq) %>% as.data.frame()
rownames(df) <- df$Var2
df$Var2 <- NULL
df %>% Heatmap()


fibs_amp2_RA$scArches_predicted <- meta_amp2_scarhces_labels_df_f$cell_type_pred

stroma_clean_h$clusters %>% unique()

fibs_amp2_RA$disease<-'amp'
stroma_clean_h$disease<-'JIA'
fibs_amp2_RA$clusters <- fibs_amp2_RA$scArches_predicted

test <- merge(fibs_amp2_RA, stroma_clean_h)

#run permutation test
test <- sc_utils(test)
prop.test <- permutation_test(test, cluster_identity = "clusters", sample_1="JIA", sample_2="amp", sample_identity="disease", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = T)
slotNames(prop.test)
results <- as.data.frame(slot(prop.test, name="results"))

