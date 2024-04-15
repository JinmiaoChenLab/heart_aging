library(Seurat)
library(pbmcapply)


##### Load results of FastIntegration and scVI #####
seurat.data = readRDS("seurat.rds")
scvi.data = read.csv("scvi.csv", row.names = 1)
colnames(scvi.data) = paste0("PCA_", 1:30)


##### Create Seurat object using raw data #####
rna.list = readRDS("rna_list.rds")
idx = split(sample(1:length(rna.list), size = length(rna.list)), 
            cut(1:length(rna.list), round(length(rna.list)/10), labels = FALSE))
rna.bind = pbmcapply::pbmclapply(
  idx, function(i) {
    rna.bind = rna.list[[i[1]]]@assays$RNA@counts
    for (j in i[2:length(i)]) {
      rna.bind = cbind(rna.bind, rna.list[[j]]@assays$RNA@counts)
    }
    return(rna.bind)
  }, mc.cores = length(idx)
)
rna.data = do.call(cbind, rna.bind)
rna.data = CreateSeuratObject(rna.data)
rna.data = NormalizeData(rna.data)

##### Append embedings of scvi and seurat into Seurat object #####
rna.data@reductions$scvi = CreateDimReducObject(as.matrix(scvi.data[Cells(rna.data),1:30]))
rna.data@reductions$seurat = CreateDimReducObject(seurat.data[Cells(rna.data),1:30])

##### Run WNN #####
rna.data = FindMultiModalNeighbors(rna.data, reduction.list = list("scvi", "seurat"), 
                                   dims.list = list(1:30, 1:30), k.nn = 30)

##### Run umap and clustering using graph of wnn #####
rna.data = RunUMAP(rna.data, nn.name = "weighted.nn", 
                   reduction.name = "umap_wnn")
rna.data = FindClusters(rna.data, graph.name = "wknn", algorithm = 2)
saveRDS(rna.data, "integrated_wnn.rds", compress = F)
