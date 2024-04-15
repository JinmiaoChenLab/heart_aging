library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(MASS)

rna.data = readRDS("snRNA.rds")
rna.data$sample = str_match(Cells(rna.data), "--(.*$)")[,2]

metadata = read.csv("metadata.txt", sep = "\t")
rownames(metadata) = metadata$Sample

metadata = metadata[rna.data$sample,]
rownames(metadata) = Cells(rna.data)
rna.data = AddMetaData(rna.data, metadata = metadata)

rna.sub = subset(rna.data, case.control == "control")
rna.sub = subset(rna.sub, Age_range_normalized != "0-9")
rna.sub = subset(rna.sub, Age_range_normalized != "10-19")
meta.sub = rna.sub@meta.data
meta.sub = meta.sub[!duplicated(meta.sub$sample),]
rownames(meta.sub) = meta.sub$sample


res = pbmclapply(
  unique(rna.sub$cell.type.group), function(ct) {
    rna = subset(rna.sub, cell.type.group == ct)
    rna$group = "young"
    rna$group[which(rna$Age_range_normalized %in% c("60-69","70-79"))] = "old"
    rna = SetIdent(rna, value = "group")
    mm = FindMarkers(rna, ident.1 = "old", ident.2 = "young", test.use = "MAST", 
                     latent.vars = c("Project", "Region", "Sex"), logfc.threshold = 0.5)
    mm$gene = rownames(mm)
    mm$cell.type = ct
    mm = mm[which(mm$p_val_adj < 0.05),]
    return(mm)
  }, mc.cores = 9
)

res = do.call(rbind, res)
res = data.frame(res)


table(res$cell.type)
saveRDS(res, "deg_two_group.rds")


a = data.frame(table(res$gene))

res.sub = res[which(res$gene %in% as.character(a$Var1[which(a$Freq > 3)])),]

plot.data = dcast(res.sub, cell.type ~ gene, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data[,1]
plot.data = plot.data[,-1]

mat_breaks <- seq(-2.5, 2.5, length.out = 100)
pheatmap::pheatmap(t(plot.data), color = colorRampPalette(c("blue", "white", "red"))(100), breaks = mat_breaks)



