library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(SeuratData)
library(FastIntegration)
library(SeuratDisk)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# metadata.txt can be downloaded from our Zenodo repository
metadata = read.csv("metadata.txt", sep = "\t")

# list of seurat objects of all samples
rna.list = readRDS("rna_list.rds")

rna.info = pbmclapply(
  rna.list, function(rna) {
    return(c(median(rna$nCount_RNA), median(rna$nFeature_RNA), str_match(Cells(rna),"--(.*)")[1,2]))
  }, mc.cores = 20
)

rna.info = do.call(rbind, rna.info)
rna.info = data.frame(rna.info)
rownames(rna.info) = rna.info$X3
colnames(rna.info) = c("ncount", "nfeature", "sample")

metadata = cbind(metadata,rna.info[metadata$Sample, 1:2])
metadata = metadata[,c("Project", "case.control", "Age_range_normalized", "Sex", "Clinical.diagnosis",
                       "Region", "source", "nfeature", "ncount")]

colnames(metadata) = c("Project", "Sample_type", "Age_range", "Sex", "Disease",
                       "Region", "Source", "nfeature", "ncount")
metadata = metadata[order(metadata$Source, metadata$Project, metadata$Age_range),]

data = metadata[,c("nfeature", "ncount")]
data$nfeature = as.numeric(data$nfeature)
data$ncount = as.numeric(data$ncount)

col.anno = metadata[,1:7]


region.color = brewer.pal(n = 6, name = 'Paired')
names(region.color) = c("LA", "RA", "LV", "RV", "AP", "S")


disease.color = brewer.pal(n = 7, name = 'Paired')
names(disease.color) = c("Normal", "ACM", "CAV", "DCM", "HCM", "MI", "Other Diagnoses")

age.color = rev(brewer.pal(n = 8, name = "RdYlBu"))
names(age.color) = sort(unique(metadata$Age_range))

project.color = brewer.pal(n = 10, name = 'Paired')
names(project.color) = sort(unique(metadata$Project))

ann_colors = list(
  Source = c(Nuclei="#2196f3", Cells="#ff5722"),
  Region = region.color,
  Disease = disease.color,
  Sex = c(F="#2196f3", M="#ff5722"),
  Age_range = age.color,
  Sample_type = c(control="#2196f3", case="#ff5722"),
  Project = project.color
)


paletteLength <- 50
myColor <- colorRampPalette(c("#2196f3", "white", "#ff5722"))(paletteLength)
myBreaks <- c(seq(500, 3000, length.out=ceiling(paletteLength/2) + 1), 
              seq(3001, 5500, length.out=floor(paletteLength/2)))

pheatmap(t(data), show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F,
         annotation_col = col.anno, annotation_colors = ann_colors, breaks = myBreaks, color = myColor)

