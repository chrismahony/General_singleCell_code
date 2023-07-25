options(bitmapType='cairo')

#############################
#read in segmented matrices, normalise and merge all sample together
setwd("/rds/projects/c/croftap-celldive01/detectionresults/matricies")

results <- dir("./", pattern = "*txt", 
    full.names = TRUE)


data = list()
for (i in 1:length(results)) {
    data[[i]] <- read.delim(results[[i]])
      
}
names(data) <- sub('.//', '', results)

library(Seurat)
library(dplyr)
library(sctransform)

#take only mean expression values for cell (remove neucle staining)
for (i in 1:length(data)) {
data[[i]] <- data[[i]] %>% select(contains('.mean'))
data[[i]] <- data[[i]] %>% select(contains('cell'))
rownames(data[[i]])<-seq_along(data[[i]][,1])
data[[i]]<-as.data.frame(t(data[[i]]))
data[[i]]<-CreateSeuratObject(data[[i]])
data[[i]]<- subset(data[[i]], subset = nCount_RNA > 0)
data[[i]]<-SCTransform(data[[i]])
}
  
#extract data
BX020_PI_EarlyRA.txt<-data[["BX020_PI_EarlyRA.txt"]]
BX028_PI_Resolver.txt<-data[["BX028_PI_Resolver.txt"]]
BX031_PI_EarlyRA.txt<-data[["BX031_PI_EarlyRA.txt"]]
BX054_PI_Resolver.txt<-data[["BX054_PI_Resolver.txt"]]
BX086_Lymphoid_EstablishedRA.txt<-data[["BX086_Lymphoid_EstablishedRA.txt"]]
BX115_Diffuse_EarlyRA.txt<-data[["BX115_Diffuse_EarlyRA.txt"]]
BX127_Diffuse_EstablishedRA.txt<-data[["BX127_Diffuse_EstablishedRA.txt"]]
BX202_Lymphoid_Resolvers.txt<-data[["BX202_Lymphoid_Resolvers.txt"]]
BX240_Lymphoid_EstablishedRA.txt<-data[["BX240_Lymphoid_EstablishedRA.txt"]]
JRP115_OA.txt<-data[["JRP115_OA.txt"]]
JRP117_OA.txt<-data[["JRP117_OA.txt"]]
JRP127_OA.txt<-data[["JRP127_OA.txt"]]

#merge all
all<-merge(x=BX054_PI_Resolver.txt, y=c(BX086_Lymphoid_EstablishedRA.txt, BX115_Diffuse_EarlyRA.txt, BX127_Diffuse_EstablishedRA.txt, BX202_Lymphoid_Resolvers.txt, BX240_Lymphoid_EstablishedRA.txt, BX031_PI_EarlyRA.txt, BX028_PI_Resolver.txt, BX020_PI_EarlyRA.txt, JRP115_OA.txt,JRP117_OA.txt, JRP127_OA.txt))


all<-FindVariableFeatures(all, assay = 'RNA')
all<-ScaleData(all, assay = 'RNA')
all<-RunPCA(all, assay = 'RNA')
all<-FindNeighbors(all, dims = 1:20)
all<-FindClusters(all, resolution = c(0.5), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.3), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.7), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.1), graph.name="RNA_snn")
all<-FindClusters(all, resolution = c(0.2), graph.name="RNA_snn")

Idents(all)<-'RNA_snn_res.0.7'
DotPlot(all, features =rownames(all)) +RotatedAxis()
table(all$RNA_snn_res.0.2)
###################################################

#rename niches based on expression of markers (for adult some were ambiguous. I.e. we consistently got a a double posative fib and vascular populaiton- this had to be extrated later, see next section)
current.sample.ids<-c( "0" , "1",  "10", "11", "12","13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "3" , "4" , "5" , "6" , "7" , "8",  "9" )
new.sample.ids<-c("mac" , "fib",  "fib", "mac", "fib","LL", "fib", "CD31_COL4A1", "fib", "Bcell", "Tcell", "pericytes", "mac", "fib", "fib", "fib" , "Tcell" , "fib" , "fib" , "mac" , "fib",  "fib" )

all$named<-all@meta.data[["RNA_snn_res.0.2"]]
                   
all@meta.data[["named"]] <- plyr::mapvalues(x = all@meta.data[["named"]], from = current.sample.ids, to = new.sample.ids)

Idents(all) <- 'named'
DotPlot(all, features =rownames(all)) +RotatedAxis()

######################################################

#divide up strange cluster
Idents(all)<-'named'
vasc<-subset(all, idents = 'vasc')
VlnPlot(vasc, features = "Cell..CD31.mean", pt.size = 0)
VlnPlot(vasc, features = "Cell..COL4A1.mean", pt.size = 0)


vasc_true<-subset(all, idents = 'vasc', subset = Cell..CD31.mean > 7 )
vasc_f<-subset(all, idents = 'vasc', subset = Cell..CD31.mean <= 7 )
ncol(vasc_true)
ncol(vasc_f)
ncol(vasc)

vasc_true<-vasc_true@meta.data
vasc_true<-subset(vasc_true, select=c("named"))
vasc_true$named<-'vasc_true'

vasc_f<-vasc_f@meta.data
vasc_f<-subset(vasc_f, select=c("named"))
vasc_f$named<-'COL4A1'

vasc_meta<-rbind(vasc_true, vasc_f)

all_meta<-all@meta.data
all_meta<-subset(all_meta, select=c("orig.ident", "named"))
all_meta<-all_meta[!rownames(all_meta) %in% rownames(vasc_meta),]
all_meta$orig.ident<-NULL
all_meta<-rbind(all_meta, vasc_meta)

nrow(all_meta)
ncol(all)

colnames(all_meta)<-"named_new"
all<-AddMetaData(all, all_meta)

##########################################################

#read in and process data coords to plot in R

data_x_y = list()
for (i in 1:length(results)) {
    data_x_y[[i]] <- read.delim(results[[i]])
    rownames(data_x_y[[i]])<-seq_along(data_x_y[[i]][,1])
    data_x_y[[i]] <- data_x_y[[i]] %>% select(c('Centroid.X.µm', 'Centroid.Y.µm'))
      }
names(data_x_y) <- sub('.//', '', results)

###########################################################

#plot one section and zoom in

#best section#
BX240_Lymphoid_EstablishedRA.txt_x_y<-data_x_y[["BX240_Lymphoid_EstablishedRA.txt"]]

Idents(all)<-'samples'
BX240<-subset(all, idents = 'BX240_Lymphoid_EstablishedRA.txt')

BX240_meta<-BX240@meta.data
BX240_Lymphoid_EstablishedRA.txt_x_y$named<-BX240_meta$RNA_snn_res.0.7

cols <- ArchR::paletteDiscrete(BX240@meta.data[, "RNA_snn_res.0.7"])
ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=0.1)+theme_classic()+scale_color_manual(values = c(cols))


BX240_Lymphoid_EstablishedRA.txt_x_y_zoom<-BX240_Lymphoid_EstablishedRA.txt_x_y[BX240_Lymphoid_EstablishedRA.txt_x_y$Centroid.X.µm < 2500,]
BX240_Lymphoid_EstablishedRA.txt_x_y_zoom<-BX240_Lymphoid_EstablishedRA.txt_x_y_zoom[BX240_Lymphoid_EstablishedRA.txt_x_y_zoom$Centroid.Y.µm > 3500,]

ggplot(BX240_Lymphoid_EstablishedRA.txt_x_y_zoom, aes(x = Centroid.X.µm, y = Centroid.Y.µm, color=named)) +
    geom_point(size=1.3)+theme_classic()+scale_color_manual(values = c(cols))

