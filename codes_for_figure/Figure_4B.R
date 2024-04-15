library(Seurat)
library(pbmcapply)
library(DISCOtoolkit)
library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(MASS)
library(sccomp)

rna.data = readRDS("snRNA.rds")

rna.sub = subset(rna.data, case.control == "control")
rna.sub = subset(rna.sub, Age_range_normalized != "0-9")
rna.sub = subset(rna.sub, Age_range_normalized != "10-19")
rna.sub$group = "young"
rna.sub$group[which(rna.sub$Age_range_normalized %in% c("60-69","70-79"))] = "old"
help(sccomp_glm)

rna.sub$Region
res = rna.sub |>
  sccomp_glm( 
    formula_composition = ~ Project + group,
    .sample =  sample, 
    .cell_group = cell.type.group, 
    bimodal_mean_variability_association = TRUE,
    cores = 10
  )
plots = plot_summary(res) 

plots$boxplot
a = plots$credible_intervals_1D
p = a[[7]]
