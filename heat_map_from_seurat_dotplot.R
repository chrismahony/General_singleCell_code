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
