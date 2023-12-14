DEGs$cluster <- "NO"
DEGs$cluster[DEGs$`Fold Change` > 0.001] <- "UP"
DEGs$cluster[DEGs$`Fold Change` < -0.001] <- "DOWN"
