library(FastIntegration)

# rna_list.rds is the list of Seurat object
# Data can be downloaded from our Zenodo repository
rna.list = readRDS("rna_list.rds") 

####### prepare data ####### 
BuildIntegrationFile(
  rna.list = rna.list, 
  tmp.dir = "/output/", # output folder
  nCores = 100, 
  nfeatures = 5000
)

####### Find anchors between samples ####### 
FastFindAnchors(tmp.dir = "/output/", nCores = 100)

####### Batch correction #######
genes = readRDS("/output/FastIntegrationTmp/raw/1.rds")
genes = rownames(genes)
idx = split(1:length(genes), cut(1:length(genes), 50, labels = FALSE))
pbmclapply(
  1:50, function(i) {
    rna.integrated = FastIntegration(tmp.dir = "/output/", npcs = 1:30, slot = "data",
                                     features.to.integrate = genes[idx[[i]]])
    saveRDS(rna.integrated, paste0("/output/FastIntegrationTmp/inte/inte_", i, ".rds"), compress = F)
  }, mc.cores = 50 
)


####### Downstream analysis #######
rna.list = pbmclapply(
  1:50, function(i){
    rna = readRDS(paste0("/output/FastIntegrationTmp/inte/inte_", i, ".rds"))
    return(rna)
  }, mc.cores = 50
)
rna.list = do.call(rbind, rna.list)
rna.list = CreateSeuratObject(rna.list)
rna.list = ScaleData(rna.list, features = rownames(rna.list))
rna.list = RunPCA(rna.list, npcs = 30, features = rownames(rna.list))
rna.list@assays$RNA@scale.data = matrix(0)


