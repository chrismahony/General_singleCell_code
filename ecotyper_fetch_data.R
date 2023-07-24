data<-stroma_clean_h@assays[["RNA"]]@data
data<-as.data.frame(data)

meta.data<-stroma_clean_h@meta.data
meta.data<-subset(meta.data, select=c(orig.ident, global_new, clusters, TYPE))

meta.data$ID<-rownames(meta.data)
meta.data$ID <- sub('-','.',meta.data$ID)
colnames(data)<-meta.data$ID

colnames(meta.data)[1] <- "Sample"
colnames(meta.data)[2] <- "CellType"

meta.data <- tibble::rownames_to_column(meta.data, "Index")
data <- tibble::rownames_to_column(data, "Gene")

library(readr)
write_delim(meta.data, "/rds/projects/c/croftap-mapjagdata/MAPJAGv2/Chris/ecotyper/stromal/annotation.txt", col_names = T, delim = "\t")
write_delim(data, "/rds/projects/c/croftap-mapjagdata/MAPJAGv2/Chris/ecotyper/stromal/data.txt",  col_names = T, delim = "\t")
