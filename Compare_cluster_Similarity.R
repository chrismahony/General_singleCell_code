library(ClusterFoldSimilarity)

Idents(obj1) <- 'my_ident'
Idents(obj2) <- 'second_ident'

list <- c(obj1, obj2)
similarity.table <- ClusterFoldSimilarity::clusterFoldSimilarity(scList = list, 
                                          sampleNames = c(1,2),
                                          topN = 1, 
                                          nSubsampling = 24)


plotClustersGraph(similarity.table)

cell.communities <- findCommunitiesSimmilarity(similarityTable = similarity.table)


similarity.table.all.values <- clusterFoldSimilarity(scList = c(obj1, obj2), 
                                                     sampleNames = c(1,2), 
                                                     topN = Inf)

similarity.table.all.values_MNP <- similarity.table.all.values
