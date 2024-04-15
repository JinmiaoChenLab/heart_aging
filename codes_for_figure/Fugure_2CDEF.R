library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(SeuratData)
library(FastIntegration)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

##### load seurat object #####
sc.rna = readRDS("scRNA.rds")
sn.rna = readRDS("snRNA.rds")

ct.shared = intersect(sc.rna$cell.type.group, sn.rna$cell.type.group)

markers = pbmclapply(
  ct.shared, function(ct) {
    sc.rna.sub = subset(sc.rna, cell.type.group == ct)
    sc.rna.sub$group = "sc"
    sn.rna.sub = subset(sn.rna, cell.type.group == ct)
    sn.rna.sub$group = "sn"
    
    rna = merge(sc.rna.sub, sn.rna.sub)
    rna = SetIdent(rna, value="group")
    rna = subset(rna, downsample = 5000)
    mm = FindMarkers(rna, ident.1 = "sn", ident.2 = "sc", group.by = "group", logfc.threshold = 1)
    mm$ct = ct
    mm$gene = rownames(mm)
    return(mm)
  }, mc.cores = 8
)

markers = do.call(rbind, markers)
markers = markers[which(markers$p_val_adj < 0.01),]

################################################################################
mm = markers
combinations = expand.grid(1:8, 1:8)
combinations$Var1 = ct.shared[combinations$Var1]
combinations$Var2 = ct.shared[combinations$Var2]

combinations$val = 0
for (i in 1:nrow(combinations)) {
  combinations$val[i] = length(intersect(
    mm$gene[which(mm$ct == combinations$Var1[i])], mm$gene[which(mm$ct == combinations$Var2[i])]
  ))
}

plot.data = dcast(combinations, Var1 ~ Var2, fill = "val")
rownames(plot.data) = plot.data$Var1
plot.data = plot.data[,-1]

for (i in 1:8) {
  plot.data[,i] = as.numeric(plot.data[,i])
}

pheatmap::pheatmap(plot.data, scale = "row")
################################################################################
mm = markers
mm$val = 1
plot.data = dcast(mm, gene ~ ct, value.var = "val", fill = 0)

pdf("deg_sc_sn_upset.pdf", width = 10, height = 6)
upset(plot.data, sets = colnames(plot.data)[2:9], 
      mb.ratio = c(0.4, 0.6),
      show.numbers = F, 
      order.by = "freq")
dev.off()


################################################################################
mm = markers
mm$avg_log2FC[which(mm$avg_log2FC > 0)] = 1
mm$avg_log2FC[which(mm$avg_log2FC < 0)] = -1

plot.data = dcast(mm, gene ~ ct, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data$gene
plot.data = plot.data[,-1]


# 235 common postive
# 263 common negative 

library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)


go_enrich <- enrichGO(gene = names(which(rowSums(plot.data) < -5)),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)



pdf("sn__neg_deg_go.pdf", width = 12, height = 6)
treeplot(pairwise_termsim(go_enrich))
dev.off()

go_enrich = go_enrich@result
go_enrich = go_enrich[which(go_enrich$p.adjust < 0.05),]
write.table(go_enrich, "sn_neg_deg_go.csv", row.names = F, col.names = T, sep = "\t", quote = F)


go_enrich <- enrichGO(gene = names(which(rowSums(plot.data) > 5)),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)
treeplot(pairwise_termsim(go_enrich))
################################################################################


mm = read.csv("sc_sn_deg.txt", sep = "\t")
mm$avg_log2FC[which(mm$avg_log2FC > 0)] = 1
mm$avg_log2FC[which(mm$avg_log2FC < 0)] = -1

plot.data = dcast(mm, gene ~ ct, value.var = "avg_log2FC", fill = 0)
rownames(plot.data) = plot.data$gene
plot.data = plot.data[,-1]

gene = names(which(rowSums(plot.data) > 5))


anno = read.csv("genes.gtf", 
                sep = "\t", skip = 5, header = F)

anno = anno[which(anno$V3 == "gene"),]


rownames(anno) = make.unique(str_match(anno$V9, "gene_name (.*?);")[,2])
anno$type = str_match(anno$V9, "gene_biotype (.*?)$")[,2]
anno$type[which(anno$type != "protein_coding")] = "ncRNA"

table(anno[gene, "type"])

table(anno[, "type"])


pie(table(anno[, "type"]), col = c('#3f51b5', '#cf2e2e'))

fisher.test(
  matrix(c(123,13503,112,19800), nrow = 2)
)

################################################################################
anno = read.csv("genes.gtf", 
                sep = "\t", skip = 5, header = F)
anno = anno[which(anno$V3 == "gene"),]
rownames(anno) = make.unique(str_match(anno$V9, "gene_name (.*?);")[,2])
anno$type = str_match(anno$V9, "gene_biotype (.*?)$")[,2]
anno$type[which(anno$type != "protein_coding")] = "ncRNA"
table(anno$type)

ct.shared = intersect(sc.rna$cell.type.group, sn.rna$cell.type.group)

ncSc = pbmclapply(
  ct.shared, function(ct) {
    rna.sub = subset(sc.rna, cell.type.group == ct)
    exp.count = rowSums(x = GetAssayData(object = rna.sub, slot = "counts") > 0)
    exp.count = exp.count/ncol(rna.sub)
    genes = names(which(exp.count > 0.01))
    return(c(table(anno[genes,"type"])["ncRNA"],table(anno[genes,"type"])["protein_coding"], ct))
  }, mc.cores = 8
)
ncSc = do.call(rbind, ncSc)
ncSc = data.frame(ncSc)
ncSc$type = "scRNA"
ncSc = data.frame(melt(ncSc))

ncSn = pbmclapply(
  ct.shared, function(ct) {
    rna.sub = subset(sn.rna, cell.type.group == ct)
    exp.count = rowSums(x = GetAssayData(object = rna.sub, slot = "counts") > 0)
    exp.count = exp.count/ncol(rna.sub)
    genes = names(which(exp.count > 0.01))
    return(c(table(anno[genes,"type"])["ncRNA"],table(anno[genes,"type"])["protein_coding"], ct))
  }, mc.cores = 8
)
ncSn = do.call(rbind, ncSn)
ncSn = data.frame(ncSn)
ncSn$type = "snRNA"
ncSn = data.frame(melt(ncSn))


plot.data = rbind(ncSc, ncSn)
plot.data = data.frame(melt(plot.data, id.vars = c("V3", "type")))
plot.data$value = as.numeric(plot.data$value)


p = ggplot(plot.data, aes(fill=variable, y=value, x=type)) + 
  geom_bar(position="fill", stat="identity") + facet_wrap(~ V3, ncol = 4)+ 
  scale_fill_manual(values=c('#cf2e2e', '#3f51b5'))
