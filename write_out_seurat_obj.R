library(DropletUtils)
#folder must not already exist i.e. you are creating a new one
write10xCounts(path="/path/", x=obj@assays$RNA@counts)
write.table(obj@meta.data %>% as.data.frame(), "/path/meta.tsv", sep="\t")
