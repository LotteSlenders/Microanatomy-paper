## 2019-11-25

# Project SCS 18 patients Experiment Rwanda

# Endothelial cells and SMC zoom

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Rwanda")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(org.Hs.eg.db)   # v3.7.0
library(tidyverse, quietly = T) # v1.2.1
library(RColorBrewer)   # v1.1-2
library(pheatmap)       # v1.0.12
library(readxl)         # 1.3.0
library(xlsx)           # 0.6.1
library(ReactomePA)     # 1.26.0

# object
load(file = "all.seur.combined.rdata")
#seuset <- all.seur.combined

#-----------------------------------------------------------------------------------------
#Extract endo cell clusters

E_clusters <- c("CD34+ Endothelial Cells", "CD34+ Endothelial Cells II")
E_cells <- WhichCells(object = all.seur.combined, ident = E_clusters)

all.seur.combined.E_clusters <- SubsetData(object = all.seur.combined, ident.use = E_clusters)

# t-SNE and Clustering
all.seur.combined.E_clusters <- RunTSNE(all.seur.combined.E_clusters, reduction.use = "pca.aligned", dims.use = 1:15, 
                                        do.fast = T, perplexity = 17)
all.seur.combined.E_clusters <- FindClusters(all.seur.combined.E_clusters, reduction.type = "pca.aligned", 
                                             resolution = 0.8, dims.use = 1:15, force.recalc = T)

TSNEPlot(all.seur.combined.E_clusters, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                           text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                           axis.text        = element_text(size = 14, face = "plain")
)


#Define marker genes per cluster
all.all.seur.combined.E_clusters.markers <- FindAllMarkers(object = all.seur.combined.E_clusters, only.pos = TRUE, min.pct = 0.25, 
                                                           thresh.use = 0.25)

#Show top2 marker genes
#all.all.seur.combined.E_clusters.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)

all.all.seur.combined.E_clusters.markers <- subset(all.all.seur.combined.E_clusters.markers, p_val_adj < 0.05)

#Save the top markers per cluster
sep.markers <- all.all.seur.combined.E_clusters.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)

#And extract the gene names
c0.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==0),"gene"]))
c1.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==1),"gene"]))
c2.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==2),"gene"]))
c3.sep.markers <- as.vector(unlist(sep.markers[which(sep.markers$cluster==2),"gene"]))

# change ident
all.seur.combined.E_clusters <- RenameIdent(object = all.seur.combined.E_clusters, old.ident.name = 0, new.ident.name = "E.0")
all.seur.combined.E_clusters <- RenameIdent(object = all.seur.combined.E_clusters, old.ident.name = 1, new.ident.name = "E.1")
all.seur.combined.E_clusters <- RenameIdent(object = all.seur.combined.E_clusters, old.ident.name = 2, new.ident.name = "E.2")
all.seur.combined.E_clusters <- RenameIdent(object = all.seur.combined.E_clusters, old.ident.name = 3, new.ident.name = "E.3")

as.character(unique(all.seur.combined.E_clusters@ident))

all.seur.combined.E_clusters <- AddMetaData(all.seur.combined.E_clusters, all.seur.combined.E_clusters@ident, "given.ident")

E_idents <- data.frame(ident = all.seur.combined.E_clusters@ident)

#-----------------------------------------------------------------------------------------
#Extract SMC clusters
#Check cells in 'bad' clusters
SMC_clusters <- c("MYH11+ Smooth Muscle Cells")
SMC_cells <- WhichCells(object = all.seur.combined, ident = SMC_clusters)

all.seur.combined.SMC_clusters <- SubsetData(object = all.seur.combined, ident.use = SMC_clusters)

# t-SNE and Clustering
all.seur.combined.SMC_clusters <- RunTSNE(all.seur.combined.SMC_clusters, reduction.use = "pca.aligned", dims.use = 1:15, 
                                          do.fast = T)
all.seur.combined.SMC_clusters <- FindClusters(all.seur.combined.SMC_clusters, reduction.type = "pca.aligned", 
                                               resolution = 0.4, dims.use = 1:15, force.recalc = T)

TSNEPlot(all.seur.combined.SMC_clusters, do.label = T, pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                             text             = element_text(family = "Arial", size = 16, face = "bold"),
                                                                                             axis.text        = element_text(size = 14, face = "plain")
)

all.seur.combined.SMC_clusters <- RenameIdent(object = all.seur.combined.SMC_clusters, old.ident.name = 0, new.ident.name = "S.0")
all.seur.combined.SMC_clusters <- RenameIdent(object = all.seur.combined.SMC_clusters, old.ident.name = 1, new.ident.name = "S.1")


as.character(unique(all.seur.combined.SMC_clusters@ident))

all.seur.combined.SMC_clusters <- AddMetaData(all.seur.combined.SMC_clusters, all.seur.combined.SMC_clusters@ident, "given.ident")

S_idents <- data.frame(ident = all.seur.combined.SMC_clusters@ident)

#-----------------------------------------------------------------------------------------
### recluster
seubset <- SubsetData(object = all.seur.combined, cells.use = c(E_cells, SMC_cells))

seubset <- RunTSNE(seubset, reduction.use = "pca.aligned", dims.use = 1:15, 
                                          do.fast = T)
seubset <- FindClusters(seubset, reduction.type = "pca.aligned", 
                                               resolution = 0.6, dims.use = 1:15, force.recalc = T)

TSNEPlot(seubset, do.label = T, pt.size = 2, label.size = 10)

given.idents <- rbind(E_idents, S_idents)
given.idents <- given.idents[rownames(as.matrix(seubset@ident)),,drop = F]
#seubset <- AddMetaData(seubset, as.character(given.idents$ident), "given.ident")
#new.ident <- as.matrix(seubset@ident)
seubset <- SetIdent(seubset, ident.use = given.idents$ident)
seubset <- AddMetaData(seubset, metadata = seubset@ident, col.name = "given.ident.2")
seubset <- SetIdent(seubset, ident.use = seubset@meta.data$res.0.6)

TSNEPlot(seubset, do.label = T, pt.size = 2, label.size = 10)

head(seubset@meta.data)

table(seubset@meta.data$res.0.6, seubset@meta.data$given.ident.2)
table(seubset@meta.data$res.0.6, seubset@meta.data$Patient)

FeaturePlot(seubset, features.plot = "MYH11")
FeaturePlot(seubset, features.plot = "CD34")
FeaturePlot(seubset, features.plot = "ACTA2")

TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "given.ident.2", plot.title = "title")

#-----------------------------------------------------------------------------------------
### plot shit
smc.markers <- readClipboard()
smc.markers <- smc.markers[smc.markers %in% rownames(seuset@data)]
smc.markers

ec.markers <- readClipboard()
ec.markers <- ec.markers[ec.markers %in% rownames(seuset@data)]
ec.markers

seubset2 <- SubsetData(object = all.seur.combined, cells.use = c(E_cells, SMC_cells))

pdf(paste(result.folder, "/", Sys.Date(), "smc.markers.pdf", sep = ""))
TSNEPlot(seubset2, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "title")
for (i in smc.markers) {
  FeaturePlot(seubset2, features.plot = i)
}
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), " ec.markers.pdf", sep = ""))
TSNEPlot(seubset2, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "title")
for (i in ec.markers) {
  FeaturePlot(seubset2, features.plot = i)
}
dev.off()

#-----------------------------------------------------------------------------------------
### DEG
seubset <- SetIdent(seubset, ident.use = given.idents$ident)

E3.markers <- FindMarkers(seubset, ident.1 = "E.3", ident.2 = c("S.0", "S.1"), only.pos = T)
E3.markers <- subset(E3.markers, p_val_adj < 0.05)

for (i in rownames(E3.markers)) {
  FeaturePlot(seubset, features.plot = i)
}

plot.the.features <- function(x){
  for (i in rownames(x)) {
    FeaturePlot(seubset, features.plot = i)
  }
}

E3.S0 <- FindMarkers(seubset, ident.1 = "E.3", ident.2 = c("S.0"), only.pos = T)
E3.S0 <- subset(E3.S0, p_val_adj <= 0.05)
E3.S1 <- FindMarkers(seubset, ident.1 = "E.3", ident.2 = c("S.1"), only.pos = T)
E3.S1 <- subset(E3.S1, p_val_adj <= 0.05)

plot.the.features(E3.S0)

E1.E2.S0 <- FindMarkers(seubset, ident.1 = c("E.1", "E.2"), ident.2 = c("S.0"), only.pos = T)
E1.E2.S0 <- subset(E1.E2.S0, p_val_adj <= 0.05)
E1.E2.S1 <- FindMarkers(seubset, ident.1 = c("E.1", "E.2"), ident.2 = c("S.1"), only.pos = T)
E1.E2.S1 <- subset(E1.E2.S1, p_val_adj <= 0.05)

E1.E2.S0.S1 <- unique(rownames(E1.E2.S0), rownames(E1.E2.S1))
rownames(E3.markers) %in% E1.E2.S0.S1

TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "nGene", plot.title = "title")
FeaturePlot(seubset, features.plot = "nGene")

plot(seubset@meta.data[,c("given.ident.2", "nGene")])
