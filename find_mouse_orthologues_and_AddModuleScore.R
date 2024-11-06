library(babelgene)

med_fibs_markers_f <- med_fibs_markers %>% filter(cluster == "FBLN1+ C5") %>% slice_head(n = 10)

genes <- orthologs(genes = med_fibs_markers_f$gene, species = "mouse", human = T)

obj <- AddModuleScore(obj, features= list(genes$symbol), name= "yourname")
