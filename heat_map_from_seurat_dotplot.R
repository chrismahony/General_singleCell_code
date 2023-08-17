library(tidyverse)

genes=c("CILP", "KIAA0040", "FOXS1", "TSPAN13", "CIITA", "CCL13", "APOL4", "CXCL9", "GBP4", "IL18BP", "IL32", "CCL5", "BIRC3", "IL4I1", "CXCL10", "MTG1", "CXCL1", "C15orf48", "CXCL6", "APOL4", "CIITA", "GBP4", "CCL8", "CXCL10", "AMTN", "GJA5", "IDO1", "GBP6", "MEOX1", "CXCL8", "CLDN1", "IL6", "CXCL2", "CXCL8", "PTGS2", "IL1B", "GBP6", "MIR3945HG", "IDO1")

Idents(merging)<-'Annotations_new_new'

dotplot<-DotPlot(merging, features = unique(genes), idents = c("periAggr_fibro", "periVasc_fibro"))

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot)

###improved colors and annotations
cols2$col
Idents(slide2_merged_f)<-'RNA_snn_res.0.3'
ha2<-columnAnnotation(show_legend = FALSE,
celltypes = levels(slide2_merged_f), 
col = list(celltypes= c("0"="#D51F26", "1"="#272E6A",  "2"="#208A42", "3"="#89288F", "4"="#F47D2B", "5"="#FEE500", "6"="#8A9FD1", "7"="#C06CAB", "8"="#D8A767")))

Heatmap(dotplot,cluster_rows = F, cluster_columns = F,top_annotation = ha2, col = circlize::colorRamp2(c(0, 0.1, 3), c('white', 'white', muted("blue"))), show_column_names = T, column_names_side = c("top"), column_names_rot = 0


