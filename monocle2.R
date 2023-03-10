
#pre process (if required)
fibs <- fibs %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(verbose = FALSE)

#choose good cluster res and caclulate markers
Idents(fibs)<-'endo_sub.cluster.0.1'

all_markers_fibs<-FindAllMarkers(fibs, only.pos = T)

all_markers_fibs %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC) -> all_markers_fibs_top50

#monocle2
library(monocle)
data <- GetAssayData(fibs[["RNA"]], slot="data")

pd <- new('AnnotatedDataFrame', data = fibs@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

data <- data[, rownames(pd)]

rna_fibros <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())


#split into groups/condition (if required)
rna_fibros_g1 <- rna_fibros[,grepl("group1", rna_fibros@phenoData@data[["PCA_groups"]], ignore.case=TRUE)]
rna_fibros_g2 <- rna_fibros[,grepl("group2", rna_fibros@phenoData@data[["PCA_groups"]], ignore.case=TRUE)]

#pseuodtime ordering and plotting
rna_fibros_g1  <- estimateSizeFactors(rna_fibros_g1)
rna_fibros_g1 <- setOrderingFilter(rna_fibros_g1, all_markers_fibs_top50$gene)
rna_fibros_g1 <- reduceDimension(rna_fibros_g1, max_components = 2, method = 'DDRTree')
rna_fibros_g1 <- orderCells(rna_fibros_g1)
rna_fibros_g1_rev <- orderCells(rna_fibros_g1, reverse = T)
plot_cell_trajectory(rna_fibros_g1, color_by = "endo_sub.cluster.0.1", cell_size = 1,show_branch_points = F)
plot_cell_trajectory(rna_fibros_g1, color_by = "TYPE", cell_size = 1,show_branch_points = F)
plot_cell_trajectory(rna_fibros_g1, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)


rna_fibros_g2  <- estimateSizeFactors(rna_fibros_g2)
rna_fibros_g2 <- setOrderingFilter(rna_fibros_g2, all_markers_fibs_top50$gene)
rna_fibros_g2 <- reduceDimension(rna_fibros_g2, max_components = 2, method = 'DDRTree')
rna_fibros_g2 <- orderCells(rna_fibros_g2)
rna_fibros_g2_rev <- orderCells(rna_fibros_g2, reverse = T)
plot_cell_trajectory(rna_fibros_g2, color_by = "endo_sub.cluster.0.1", cell_size = 1,show_branch_points = F)
plot_cell_trajectory(rna_fibros_g2, color_by = "TYPE", cell_size = 1,show_branch_points = F)
plot_cell_trajectory(rna_fibros_g2, color_by = "Pseudotime", cell_size = 1,show_branch_points = F)
