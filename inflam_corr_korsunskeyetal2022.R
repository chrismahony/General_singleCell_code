#create meta data splot
med_fibros$sample_tissue_cluster<-paste(med_fibros$DonorID, med_fibros$Tissue, med_fibros$Cluster_name, sep=".")


paste(med_fibros$InflamScore, med_fibros$LibraryID, sep="_") %>% table()

#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue_cluster"
dotplot<-DotPlot(med_fibros, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"




dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$LibraryID, med_fibros$InflamScore, sep="_") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep="_")


index <- match(dotplot_data$id_1, inflam_df$._1)
dotplot_data$inflam <- inflam_df$._2[index]


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(dotplot_data, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


ml = lm(Runx1~InflamScore, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)
ggscatter(dotplot_data, x = "Runx1", y = "InflamScore",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~condition_2)



#filter if you want a specific cluster
dotplot_data_C4 <- dotplot_data[dotplot_data$id_3 == "SPARC+COL3A1+ C4",]
dotplot_data_C4 <- dotplot_data_C4[dotplot_data_C4$id_2 == "SalivaryGland",]

ggplot(dotplot_data_C4, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)
