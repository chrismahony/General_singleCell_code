library(EnhancedVolcano)

genes=c("gene1", "gene2")

EnhancedVolcano(df,
    lab = rownames(df),
    x = 'avg_log2FC',
    y = 'p_val_adj',
        selectLab = genes,
    title = 'title',
    subtitle = "GEX, red=p_adj<0.05 & FC > 0.25",
    pCutoff = 0.05,
    FCcutoff = 0.25,
    pointSize = 3.0,
    labSize = 3,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
      legendPosition = 'none',
    legendLabSize = 10,
    legendIconSize =3.0,
    legendLabels=c('NS','p<0.05 & FC > 0.25'),
    col=c('black', 'black', 'black', 'red3'),
          drawConnectors = TRUE,
    widthConnectors = 0.3,
     xlab = bquote(~Log[2]~ 'fold change'),
    boxedLabels = T,
    borderWidth = 0.5
    )
