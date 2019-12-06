## 2019-6-24

# Project SCS 18 patients Experiment Oman

# Find candidates for staining SMC derived foamy MCs

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Oman")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.8
library(org.Hs.eg.db)   # v3.7.0
library(SingleR)        # v0.2.2

# object
seuset <-  readRDS(file = "Myeloid macrophage clusters.RDS")
load(file = "all.seur.combined.rdata")

# set seed
set.seed(31)

# check tSNE
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "check")

#-----------------------------------------------------------------------------------------
# check stuff
VlnPlot(seuset, features.plot = c("ACTA2", "MYH11", "CD68"))
DoHeatmap(seuset, genes.use = c("ACTA2", "MYH11", "CD68"))
DoHeatmap(seuset, genes.use = c("ACTA2", "MYH11", "CD68"), group.by = "Sex")

#-----------------------------------------------------------------------------------------
# what differentiates cluster 2 from the others?

all.markers.strict <-  FindAllMarkers(seuset, logfc.threshold = 0.6,
                              min.pct = 0.25)

all.markers.default <- FindAllMarkers(seuset)

write.table(all.markers.strict, file = paste(result.folder, "/", Sys.Date(), "all_markers_strict.txt", sep = ""))
write.table(all.markers.default, file = paste(result.folder, "/", Sys.Date(), "all_markers_default.txt", sep = ""))

markers2.strict <- subset(all.markers.strict, cluster == 2 & p_val_adj <= 0.05)
markers2.default <- subset(all.markers.default, cluster == 2 & p_val_adj <= 0.05)

DoHeatmap(seuset, genes.use = c(markers2.strict$gene[1:30], "ACTA2", "MYH11","CD68"))
DoHeatmap(seuset, genes.use = c(markers2.strict$gene[31:50], "ACTA2", "MYH11","CD68"))
DoHeatmap(seuset, genes.use = c(markers2.strict$gene[51:70], "ACTA2", "MYH11","CD68"))
DoHeatmap(seuset, genes.use = c(markers2.default$gene[1:20], "ACTA2", "MYH11"))

# 1 - 30
pdf(paste(result.folder, "/", Sys.Date(), " Strict markers 1-30.pdf", sep = ""))
for(i in markers2.strict$gene[1:30]){
FeaturePlot(all.seur.combined, features.plot = i)
}
dev.off()

# 31 - 50
pdf(paste(result.folder, "/", Sys.Date(), " Strict markers 31-50.pdf", sep = ""))
for(i in markers2.strict$gene[31:50]){
  FeaturePlot(all.seur.combined, features.plot = i)
}
dev.off()

# 51 - 70
pdf(paste(result.folder, "/", Sys.Date(), " Strict markers 51-70.pdf", sep = ""))
for(i in markers2.strict$gene[51:70]){
  FeaturePlot(all.seur.combined, features.plot = i)
}
dev.off()

VlnPlot(seuset, features.plot = c("ACTA2", "MMP9", "APOC1"))
FeaturePlot(all.seur.combined, features.plot = "CYP27A1")

# group by patient for SPP1 (calcification)
DoHeatmap(seuset, genes.use = c("SPP1", "ACTA2", "CD68"), group.by = "Patient", group.label.rot = T)
# group by sex
DoHeatmap(seuset, genes.use = c("SPP1", "ACTA2", "CD68"), group.by = "Sex", group.label.rot = T)

# possible candidates with low expression in SMCs
DoHeatmap(seuset, genes.use = c("ACTA2", "CD68","SPP1", "APOC1", "MMP9", "SCD", "LHFPL2", "IL4I1", "TREM2"))
DoHeatmap(seuset, genes.use = c("ACTA2", "CD68","SPP1", "APOC1", "MMP9", "SCD", "LHFPL2", "IL4I1", "TREM2"), group.by = "Sex", group.label.rot = T)

# possible candidates from literature
DoHeatmap(seuset, genes.use = c("VCAM1","CCL2", "ABCA1", "LGALS3", "ITGAM", "ADGRE1", "TAGLN"))


#--------------------------------
# MMPs

to.plot <- c("ACTA2", "CD68", "MMP1", "MMP2", "MMP3","MMP7", "MMP8", "MMP9", "MMP19")
pdf(paste(result.folder, "/", Sys.Date(), " MMPs.pdf", sep = ""))
for(i in to.plot){
  FeaturePlot(all.seur.combined, features.plot = i)
}
dev.off()


