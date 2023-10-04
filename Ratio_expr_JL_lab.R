#make a new meta.data column
med_fibros$cluster_condition_sample <- paste(med_fibros$Cluster_name, med_fibros$Case, med_fibros$LibraryID, sep="_")
Idents(med_fibros)<-'cluster_condition_sample'

#choose the idents that you want to plot, here I sub set clusters and then tissue/conditon
idents <- levels(med_fibros)[grep('C2_|C6_', levels(med_fibros))]
idents <- idents[grep('GutControl|GutNonisnflamed', idents)]

#dotplot and puss avg.scaled.expr
dotplot<-DotPlot(med_fibros, features = c("COL1A1", "COL3A1"), idents = idents)

dotplot<-dotplot$data

library(dplyr)
library(tidyr)
dotplot<-dotplot %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

#tidy thihgs up a bit
dotplot <- t(dotplot) %>% as.data.frame()
colnames(dotplot) <-  dotplot[1,]
dotplot <- dotplot[-1,]
dotplot$COL1A1 <- as.double(dotplot$COL1A1)
dotplot$COL3A1 <- as.double(dotplot$COL3A1)
dotplot$ratio <- as.double(dotplot$COL1A1/dotplot$COL3A1)
dotplot$meta <- rownames(dotplot)
library(splitstackshape)
dotplot <- cSplit(dotplot, splitCols = "meta", sep="_")

#plot and stats
library(ggpubr)
ggplot(dotplot, aes(x=ratio, y=factor(meta_2, levels=c("GutControl","GutNonisnflamed")))) +   geom_violin(trim=F,)+ coord_flip()+ geom_boxplot(width=0.1)+geom_jitter(shape=16, position=position_jitter(0.2))+ scale_fill_grey() + theme_classic()

compare_means(ratio ~ meta_2, data = dotplot)


#a niver box and whiskers plot
wisteria <- c("grey65", "burlywood3", "khaki2", "plum1", "lightcyan2", "cornflowerblue", "slateblue3")
 
  ggplot(dotplot, aes(x=factor(meta_2, levels=c("GutControl","GutNonisnflamed")), y=ratio)) +
  geom_boxplot(aes(fill = meta_2),            #You can make a box plot too!
               alpha = 0.8, width = 0.7) +      
  geom_point(aes(fill = meta_2), shape = 21, color = "black", alpha = 0.8,
             position = position_jitter(width = 0.1, seed = 666))+
  scale_fill_manual(values = wisteria[c(3, 1, 6)]) +
  labs(x = "Species",
       y = "Bill length") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line = element_line(size = 1.2),
        text = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold")
        )
