
#pull GEX from seurat obj (mine is all)
gex <- t(all@assays[["RNA"]]@data) %>% as.data.frame()
head(gex)

#add any meta.data you want to plot
gex$named <- all$named_niches
gex$condition <- all$condition

#get colors (optional)
cols <- as.data.frame(ArchR::paletteDiscrete(all@meta.data[, "named_niches"]))
colnames(cols) <- "colors"

library(tidyverse)

#plot. No outliers and no whiskers
gex %>% 
  pivot_longer(cols = gene, 
               names_to = "gene", 
               values_to = "expression") %>% 
ggplot(aes(x=gene, y=expression, fill=named))+
    geom_boxplot(outlier.shape = NA, coef = 0)+facet_wrap(~ condition)+theme_minimal()+ylim(0,2000)+theme(axis.ticks = element_blank())+scale_fill_manual(values=cols$colors)


#plot. With outliers and with whiskers
gex %>% 
  pivot_longer(cols = gene, 
               names_to = "gene", 
               values_to = "expression") %>% 
ggplot(aes(x=gene, y=expression, fill=named))+
    geom_boxplot()+facet_wrap(~ condition)+theme_minimal()+ylim(0,2000)+theme(axis.ticks = element_blank())+scale_fill_manual(values=cols$colors)
