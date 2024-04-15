library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(SeuratData)
library(FastIntegration)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(treemap)

##### load seurat object #####
sc.rna = readRDS("scRNA.rds")
sn.rna = readRDS("snRNA.rds")


##### The following codes are for scRNA-seq data #####
plot.data = data.frame(
  subgroup = unique(sc.rna$cell.type),
  value = as.numeric(table(sc.rna$cell.type)[unique(sc.rna$cell.type)])
)
meta = sc.rna@meta.data
meta = meta[!duplicated(meta$cell.type),]
rownames(meta) = meta$cell.type
plot.data$group = meta[plot.data$subgroup, "cell.type.group"]


cell.color = read.csv("cell_type.txt",sep = "\t") # color file, can be found in our Zenodo repository
plot.data$mycolors =  cell.color[plot.data$subgroup,1]
treemap(plot.data, index=c("group","subgroup"), vSize="value", type="color",
        fontsize.labels=c(17,12),
        fontcolor.labels=c("black","white"),    
        fontface.labels=c(2,1),
        bg.labels=c("transparent"),            
        align.labels=list(
          c("center", "center"),
          c("right", "bottom")
        ),
        overlap.labels=0.5,                     
        inflate.labels=F, 
        palette = cell.color[plot.data$subgroup,1],
        vColor ="mycolors", title = "snRNA-seq"
)
