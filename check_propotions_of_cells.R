obj$sample_condition<-paste(obj$sample, obj$condition, sep=".") #condition would be resting/peak/resolved etc

pt <- table(obj$sample_condition, obj$cluster.name)
pt <- as.data.frame(pt)

library(tidyverse)
pt_cluster1 = subset(pt, Var2 == "cluster1")
pt_cluster2 = subset(pt, Var2 == "cluster2")
pt_cluster3 = subset(pt, Var2 == "cluster3")

colnames(pt_cluster1) <- c("Var1", "Var2", "cluster1")
colnames(pt_cluster2)<- c("Var1", "Var2", "cluster2")
colnames(pt_cluster3)<- c("Var1", "Var2", "cluster3")

pt_master<-pt_cluster1
pt_master$cluster2<-cluster2$cluster2
pt_master$cluster3<-cluster3$cluster3

pt_master <- pt_master %>% select(-one_of('Var2'))
library(tidyverse)
pt_master <- pt_master %>% remove_rownames %>% column_to_rownames(var="Var1")
pt_master <- pt_master/rowSums(pt_master)
pt_master <- as.data.frame(pt_master)

pt_master$condition<-rownames(pt_master)

library(splitstackshape)
pt_master<-cSplit(pt_master, splitCols="condition", sep=".")

#chnage the levels for what you have i your object. They will apear in the order they are written
cluster1 <- ggplot(pt_master, aes(x=factor(condition, level=c('control', 'initiation', 'peak', 'resolving', 'resolved', 'persistant')), y=cluster1)) + 
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ggtitle("Chodl")+theme_classic()+RotatedAxis()
cluster1
res.aov_cluster1 <- aov(cluster1 ~ condition, data = pt_master)
summary(res.aov_cluster1)
stats_cluster1 <- TukeyHSD(res.aov_cluster1)
stats_cluster1

#repeat for cluster 2,3 etc

