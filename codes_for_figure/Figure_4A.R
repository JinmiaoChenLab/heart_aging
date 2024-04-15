library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)


rna.data = readRDS("snRNA.rds")

rna.sub = subset(rna.data, case.control == "control")
rna.sub = subset(rna.sub, Age_range_normalized != "0-9")
rna.sub = subset(rna.sub, Age_range_normalized != "10-19")

rna.sub$group = "young"
rna.sub$group[which(rna.sub$Age_range_normalized %in% c("60-69","70-79"))] = "old"

meta.sub = rna.sub@meta.data
meta.sub = meta.sub[!duplicated(meta.sub$sample),]
rownames(meta.sub) = meta.sub$sample

plot.data = data.frame(table(rna.sub$sample, rna.sub$cell.type.group))
names(plot.data) = c("sample", "cell_type", "count")
plot.data$sample[1:5]
plot.data$sample = as.character(plot.data$sample)
plot.data$cell_type = as.character(plot.data$cell_type)
plot.data$sample[1:5]

plot.data$ratio = 0
for (i in 1:nrow(plot.data)){
  plot.data$ratio[i] = plot.data$count[i]/sum(plot.data$count[which(plot.data$sample == plot.data$sample[i])])
}

plot.data$age = meta.sub[plot.data$sample, "group"]
plot.data$age = factor(plot.data$age, levels = c("young", "old"))

plot.list = list()
my_comparisons <- list( c("young", "old"))

for (i in 1:length(unique(plot.data$cell_type))) {
  plot.data.sub  = plot.data[which(plot.data$cell_type == unique(plot.data$cell_type)[i]),]
  plot.list[[i]] = ggboxplot(plot.data.sub, "age", "ratio",
                             color = "age", palette = c("#fbbc04", "#c53929"),
                             add = "jitter")  + NoLegend() + ggtitle(unique(plot.data$cell_type)[i]) +
    stat_compare_means(comparisons = my_comparisons)
}

cowplot::plot_grid(plotlist = plot.list, ncol = 3)

