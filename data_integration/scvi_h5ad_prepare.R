library(Seurat)
library(pbmcapply)
library(SeuratData)
library(stringr)
library(SeuratDisk)
options(bitmapType = "cairo")

rna.list = readRDS("rna_list.rds")
genes = SelectIntegrationFeatures(rna.list, nfeatures = 5000)


##### subset data of feature gene #####
for (i in 1:length(rna.list)) {
  rna.list[[i]] = subset(rna.list[[i]], features = genes)
}


##### Merge rna list into one Seurat Object #####
idx = split(sample(1:length(rna.list), size = length(rna.list)), 
            cut(1:length(rna.list), round(length(rna.list)/10), labels = FALSE))
rna.bind = pbmcapply::pbmclapply(
  idx, function(i) {
    rna.bind = GetAssayData(object = rna.list[[i[1]]], assay = "RNA", slot = "counts")
    for (j in i[2:length(i)]) {
      rna.bind = cbind(
        rna.bind, 
        GetAssayData(object = rna.list[[j]], assay = "RNA", slot = "counts")
      )
    }
    return(rna.bind)
  }, mc.cores = length(idx)
)
rna.bind = do.call(cbind, rna.bind)
rna.bind = CreateSeuratObject(rna.bind)


##### Convert Seurat Object to h5ad file#####
rna.bind$sample = str_match(Cells(rna.bind), "--(.*$)")[,2]
rna.bind@assays$RNA@meta.features[,"gene"] = rownames(rna.bind) 
SaveH5Seurat(rna.bind, filename = "raw_rna.h5seurat", overwrite = T)
Convert("raw_rna.h5seurat", dest = "h5ad", overwrite = T)
