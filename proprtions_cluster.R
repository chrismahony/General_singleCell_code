pt <- table(myobj$cluster, myobj$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())



##### better plot

pt2 <- table(stia_ecotyper$cluster.name, stia_ecotyper$State)
pt2 <- as.data.frame(pt2)
pt2$Var1 <- as.character(pt2$Var1)

ggplot(pt2, aes(x = Var2, y = Freq, fill = Var1))  + 
        geom_bar(stat = 'identity', position = position_fill())+
        theme(axis.text.y= element_blank(), axis.ticks.y = element_blank()) +
        coord_flip() +  scale_fill_manual(values = cols$colors)+
        theme(legend.position = 'bottom', panel.border = element_rect(colour = "black", fill=NA, size=0.5), panel.background = element_blank()) + 
        labs(fill = 'Tissue-Defined Cluster', y = 'Cluster Frequency in Sample') + 
        guides(fill = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 2)), alpha = FALSE)

