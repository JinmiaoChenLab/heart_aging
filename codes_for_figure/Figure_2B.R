library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(SeuratData)
library(FastIntegration)
library(SeuratDisk)
library(ggplot2)


##### load seurat object #####
sc.rna = readRDS("scRNA.rds")
sn.rna = readRDS("snRNA.rds")

sc.average = AverageExpression(sc.rna, group.by = "cell.type.group")
sn.average = AverageExpression(sn.rna, group.by = "cell.type.group")

sc.average = sc.average$RNA
sn.average = sn.average$RNA


##### only use the overlapped feature genes #####
feature1 = readRDS("./scRNA/FastIntegrationTmp/others/features.rds")
feature2 = readRDS("./snRNA/FastIntegrationTmp/others/features.rds")
features = intersect(feature1, feature2)


plot.data = cor(sc.average[features,], sn.average[features,], method = "spearman")
pheatmap::pheatmap(t(plot.data), display_numbers = T)



