library(scProportionTest)

#create object for test
test <- sc_utils(seurat_obj)

#check idents names (optional)
table(seurat_obj$condition)

#run test and plot. "cluster.name" is the name of my clusters (e.g. fibroblast, macrophage etc), sample_1/smaple_2 is the name of the conditions you are comparing (e.g. control treated), smaple identity would be somehting like treatments
prop.test <- permutation_test(test, cluster_identity = "cluster.name", sample_1="control", sample_2="initiation", sample_identity="condition", n_permutations=10000)
permutation_plot(prop.test, FDR_threshold = 0.01, log2FD_threshold = 0.58, order_clusters = F)
