genes <- c("gene1", "gene2", "gene3")

sobj <- AddModuleScore(obj, features = list(genes), name="module")

#alway add a 1 to the module name
FeaturePlot(obj, "module1")

#play around with max/min cutoff (q10=bottom 10%)
FeaturePlot(obj, "module1", max.cutoff = "q90", min.cutoff = "q10")


#coorelate gene module with a speific gene
obj$sample_condition<-paste(obj$sample, obj$condition, sep=".")  #condition could be tisse type?

Idents(obj)<-"sample_condition"
dotplot<-DotPlot(obj, features = "module1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "module1"

dotplot_MMP1<-DotPlot(obj, features = "Mmp1")
dotplotmmp1_data<-dotplot_Mmp1[["data"]]
dotplotmmp1_data <- subset(dotplotmmp1_data, select = c(avg.exp.scaled, id))
names(dotplotmmp1_data)[names(dotplotmmp1_data)=="avg.exp.scaled"] <- "Mmp1"

dotplot_data$Mmp1<-dotplotmmp1_data$Mmp1
dotplot_data$condition<-dotplot_data$id
library(splitstackshape)
dotplot_data<-cSplit(dotplot_data, splitCols="condition", sep=".")


ggplot(dotplot_data, aes(x = Mmp1, y = Module1)) +
    geom_point(aes(color = factor(condition_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none")


