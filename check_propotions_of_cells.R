library(tidyr)
obj$sample_condition<-paste(obj$orig.ident, obj$condition, sep=".") #condition would be resting/peak/resolved etc

pt <- table(obj$sample_condition, obj$cluster.name)
pt <- as.data.frame(pt)

pt<-pt %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  as.data.frame() 

library(tidyverse)
pt <- pt %>% remove_rownames %>% column_to_rownames(var="Var1")
pt <- pt/rowSums(pt)

pt$condition<-rownames(pt)

library(splitstackshape)
pt<-cSplit(pt, splitCols="condition", sep=".")


cluster1 <- ggplot(pt, aes(x=factor(condition_2, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=cluster1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Cluster1")+theme_classic()+RotatedAxis()
cluster1
res.aov_cluster1 <- aov(cluster1` ~ condition_2, data = pt)
summary(res.aov_cluster1)
stats_cluster1 <- TukeyHSD(res.aov_cluster1)
stats_cluster1

#repeat for cluster 2,3 etc

