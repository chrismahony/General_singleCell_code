library(DropletUtils)
write10xCounts(path="/path/", x=obj@assays$RNA@counts)
