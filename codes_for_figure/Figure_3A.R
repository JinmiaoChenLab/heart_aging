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

rna.sub = subset(rna.data, case.control == "control")
rna.sub = subset(rna.sub, Age_range_normalized != "<10")
rna.sub = subset(rna.sub, Age_range_normalized != "10-19")

rna.sub = subset(rna.sub, Project == "ERP123138")

meta.sub = rna.sub@meta.data
meta.sub = meta.sub[!duplicated(meta.sub$sample),]
rownames(meta.sub) = meta.sub$sample


pbmclapply(
  unique(rna.sub$cell.type.group), function(ct) {
    rna = subset(rna.sub, cell.type.group == ct)
    
    exp.count = rowSums(x = GetAssayData(object = rna, slot = "counts") > 0) 
    exp.count = exp.count/ncol(rna)
    genes = names(which(exp.count > 0.01))
    
    rna.average = AverageExpression(rna, group.by = "sample")
    rna.average = rna.average$RNA
    
    meta = rna@meta.data
    meta = meta[!duplicated(meta$sample),]
    rownames(meta) = meta$sample
    meta = meta[colnames(rna.average),]
    meta$Age_range_normalized = factor(meta$Age_range_normalized, levels = c("20-29", "30-39", "40-49", "60-69", "70-79"))
    
    cor.res = pbmclapply(
      genes, function(g) {
        a = lm(rna.average[g,] ~ as.numeric(meta$Age_range_normalized) + meta$Sex)
        a = summary(a)
        return(as.numeric(a$coefficients[2,]))
      }, mc.cores = 10
    )
    cor.res = do.call(rbind, cor.res)
    cor.res = data.frame(cor.res)
    colnames(cor.res) = c("est", "error", "t", "p")
    cor.res$gene = genes
    cor.res$cell.type = ct
    
    saveRDS(cor.res, paste0(sub("/", " ", ct), "_ERP123138.rds"))
  }, mc.cores = 10
)

###############################################

metadata = read.csv("metadata.txt", sep = "\t")
rownames(metadata) = metadata$Sample

metadata = metadata[rna.data$sample,]
rownames(metadata) = Cells(rna.data)
rna.data = AddMetaData(rna.data, metadata = metadata)

rna.sub = subset(rna.data, case.control == "control")
rna.sub = subset(rna.sub, Age_range_normalized != "<10")
rna.sub = subset(rna.sub, Age_range_normalized != "10-19")


deg = readRDS("deg_v2.rds")

rna = subset(rna.sub, Project == "EGAS00001006374")
deg.sub = deg[which(deg$project == "EGAS00001006374"),]
deg.sub = deg.sub[which(deg.sub$p < 0.01),]

age.score.1 = pbmclapply(
  unique(deg.sub$cell.type), function(ct) {
    a = deg.sub[which(deg.sub$cell.type == ct),]
    rna = AddModuleScore(
      object = rna,
      ctrl = 30,
      features = list(a$gene[which(a$est > 0)]),
      name = "age_positive",
      seed = 1
    )
    rna = AddModuleScore(
      object = rna,
      ctrl = 30,
      features = list(a$gene[which(a$est < 0)]),
      name = "age_negative",
      seed = 1
    )
    a1 = data.frame(pos = rna$age_positive1, age = rna$Age_range_normalized, 
                    cell.type = ct, project = "EGAS00001006374", type = "pos")
    a2 = data.frame(pos = rna$age_negative1, age = rna$Age_range_normalized, 
                    cell.type = ct, project = "EGAS00001006374", type = "neg")
    a = rbind(a1,a2)
    return(a)
  }, mc.cores = 9
)
age.score.1 = do.call(rbind, age.score.1)

rna = subset(rna.sub, Project == "ERP123138")
deg.sub = deg[which(deg$project == "ERP123138"),]
deg.sub = deg.sub[which(deg.sub$p < 0.01),]

age.score.2 = pbmclapply(
  unique(deg.sub$cell.type), function(ct) {
    a = deg.sub[which(deg.sub$cell.type == ct),]
    rna = AddModuleScore(
      object = rna,
      ctrl = 30,
      features = list(a$gene[which(a$est > 0)]),
      name = "age_positive",
      seed = 1
    )
    rna = AddModuleScore(
      object = rna,
      ctrl = 30,
      features = list(a$gene[which(a$est < 0)]),
      name = "age_negative",
      seed = 1
    )
    a1 = data.frame(pos = rna$age_positive1, age = rna$Age_range_normalized, 
                    cell.type = ct, project = "ERP123138", type = "pos")
    a2 = data.frame(pos = rna$age_negative1, age = rna$Age_range_normalized, 
                    cell.type = ct, project = "ERP123138", type = "neg")
    a = rbind(a1,a2)
    return(a)
  }, mc.cores = 9
)
age.score.2 = do.call(rbind, age.score.2)
age.score = rbind(age.score.1, age.score.2)


age.score.summary = age.score %>%
  group_by(cell.type, project, age, type) %>%
  summarise_at(vars(pos), list(name = median))

age.score.summary = data.frame(age.score.summary)
colnames(age.score.summary)[5] = "score"


age.score.summary$type = paste0(age.score.summary$project, "_", age.score.summary$type)

a = age.score.summary %>%
  group_by(type, cell.type) %>%
  summarise_at(vars(score), list(name = min))
a = data.frame(a)

b = age.score.summary %>%
  group_by(type, cell.type) %>%
  summarise_at(vars(score), list(name = max))
b = data.frame(b)

a = data.frame(min = a$name, max = b$name, row.names = paste0(b$type,b$cell.type))

age.score.summary$score[grep("pos",age.score.summary$type)] = age.score.summary$score[grep("pos",age.score.summary$type)] - a[paste0(age.score.summary$type[grep("pos",age.score.summary$type)],age.score.summary$cell.type[grep("pos",age.score.summary$type)]),"min"]

age.score.summary$score[grep("neg",age.score.summary$type)] = age.score.summary$score[grep("neg",age.score.summary$type)] - a[paste0(age.score.summary$type[grep("neg",age.score.summary$type)],age.score.summary$cell.type[grep("neg",age.score.summary$type)]),"max"]


ggplot(data=age.score.summary, aes(x=age, y=score, group=type)) +
  geom_line(aes(color=type))+
  geom_point(aes(color=type)) + facet_wrap(~ cell.type, ncol = 5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


