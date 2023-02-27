
obj$sample_condition<-paste(obj$sample, obj$condition, sep=".")

Idents(obj)<-"sample_condition"
dotplot<-DotPlot(obj, features = "Thy1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"

dotplot_MMP1<-DotPlot(obj, features = "Mmp1")
dotplotmmp1_data<-dotplot_Mmp1[["data"]]
dotplotmmp1_data <- subset(dotplotmmp1_data, select = c(avg.exp.scaled, id))
names(dotplotmmp1_data)[names(dotplotmmp1_data)=="avg.exp.scaled"] <- "Mmp1"


dotplot_data$Mmp1<-dotplotmmp1_data$Mmp1
dotplot_data$condition<-dotplot_data$id
library(splitstackshape)
dotplot_data<-cSplit(dotplot_data, splitCols="condition", sep=".")


ggplot(dotplot_data, aes(x = Mmp1, y = Thy1)) +
    geom_point(aes(color = factor(condition_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none")


ml = lm(Mmp1~Thy1, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)
ggscatter(dotplot_data, x = "Mmp1", y = "Thy1",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)
