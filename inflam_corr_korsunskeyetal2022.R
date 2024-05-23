#check inflam scoire with GEX per sample in each tissue

med_fibros$sample_tissue<-paste(med_fibros$SampleID, med_fibros$Tissue, sep=".")


#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue"
dotplot<-DotPlot(med_fibros, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


ml = lm(Runx1~InflamScore, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)
ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)




###############################
#check inflam score with GEX in each sample/cluster

#create meta data splot
med_fibros$sample_tissue_cluster<-paste(med_fibros$SampleID, med_fibros$Tissue, med_fibros$Cluster_name, sep=".")


paste(med_fibros$InflamScore, med_fibros$LibraryID, sep="_") %>% table()

#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros)<-"sample_tissue_cluster"
dotplot<-DotPlot(med_fibros, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)


ml = lm(Runx1~InflamScore, data = dotplot_data)
summary(ml)$r.squared

library(ggpubr)
ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)



#filter if you want a specific cluster
dotplot_data_C4 <- final_df[final_df$id_3 == "SPARC+COL3A1+ C4",]
dotplot_data_C4 <- dotplot_data_C4[dotplot_data_C4$id_2 == "Synovium",]

ggplot(dotplot_data_C4, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_2))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") +
        facet_wrap(~id_2)




ml = lm(Runx1~InflamScore, data = dotplot_data_C4)
summary(ml)$r.squared

library(ggpubr)

dotplot_data_syn <- final_df[final_df$id_2 == "Synovium",]

ggscatter(dotplot_data_C4, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = 0, label.y = 2)+
        facet_wrap(~id_2)


