library(ggtree)

Idents(stia2021_rna) <- "cluster.name"

rds <- BuildClusterTree(stia2021_rna, 
                        dims = 1:30)

data.tree <- Tool(object = rds, 
                  slot = "BuildClusterTree")



gg_tr <- ggtree(data.tree, 
                layout = "rectangular",
                alpha = .1,
                size = 5) +
  geom_tiplab(hjust = 1.1, size = 3) +
  geom_treescale() +
  theme(plot.margin = unit(c(1,5,1,1), "mm")) 

plot(gg_tr)
