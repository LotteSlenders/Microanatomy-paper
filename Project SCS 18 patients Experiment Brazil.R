### 2019-3-28

# Project SCS 18 patients Experiment Brazil

# make plots

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Brazil")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(org.Hs.eg.db)   # v3.7.0
library(SingleR)        # v0.2.2

# object
seuset <- readRDS("seuset 18 patients 16 communities2.RDS")
seubset <- readRDS("seuset 18 patients 16 communities2.RDS")

#-----------------------------------------------------------------------------------------
### 2019-4-2 Virginia plots Michal
## subset SMCs
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "ident")

seubset <- SubsetData(seuset,
                      ident.use = 8)
seubset <- FindVariableGenes(seubset)
seubset <- RunPCA(seubset, do.print = F)
seubset <- FindClusters(object = seubset, 
                        reduction.type = "pca", 
                        resolution = 1, 
                        dims.use = 1:10, 
                        force.recalc = T,
                        random.seed = 91 )
seubset <- RunTSNE(seubset, dims.use = 1:10)
TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "ident")

singler <- CreateSinglerObject(counts = seubset@raw.data[,rownames(seubset@meta.data)],       
                               annot = seubset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "SMCs", 
                               min.genes = 0,                          # if starting from seurat, set min.genes to 0 
                               technology = "CEL-seq2", 
                               species = "Human", 
                               citation = "",
                               ref.list = list(), 
                               normalize.gene.length = FALSE, 
                               variable.genes = "de",
                               fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
                               do.signatures = TRUE,
                               clusters = seubset@ident,
                               do.main.types = TRUE, 
                               reduce.file.size = TRUE, 
                               numCores = 4
)


pdf(paste(result.folder, "/", Sys.Date(), " Virginia genes SMCs.pdf", sep = ""))
TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "ident")
#SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 20)

# Myh11+ groups: (+ for Myh11, Acta2, Cnn1, Tagln)
to.plot <- c("MYH11", "ACTA2", "CNN1", "TAGLN", "CD200", "DES")

DoHeatmap(seubset,
          genes.use = to.plot,
          title = "MYH11+ groups")
for (i in to.plot) {
  FeaturePlot(seubset, features.plot = i)
}
lapply(to.plot, function(x){
  VlnPlot(seubset, features.plot = x)
})

# Transitional/stem like groups:
to.plot <- c("VCAM1", "CD34", "PDGFRA", "ID4", "IL6", "SPP1", "DCN", "LUM", "IL4R")

DoHeatmap(seubset,
          genes.use = to.plot,
          title = "Transitional/stem like groups")
for (i in to.plot) {
  FeaturePlot(seubset, features.plot = i)
}
lapply(to.plot, function(x){
  VlnPlot(seubset, features.plot = x)
})

# Chondrocyte-like groups: 
to.plot <- c("SOX9", "RUNX2", "TRPV4", "S100B")

DoHeatmap(seubset,
          genes.use = to.plot,
          title = "Chondrocyte-like groups")
lapply(to.plot, function(x){
  VlnPlot(seubset, features.plot = x)
})


dev.off()
#-----------------------------------------------------------------------------------------
t <- seuset@raw.data[,rownames(seuset@meta.data)]
saveRDS(t, "Alona file.RDS")

#-----------------------------------------------------------------------------------------
# 2019-4-15 fib vs athero

fib.vs.athero.endo <- read.delim("Fib vs Athero Endo not transformed.txt")
fib.vs.athero.SMC <- read.delim("Fib vs Athero SMC not transformed.txt")

# transform
fib.vs.athero.endo <- log2(fib.vs.athero.endo + 0.1)
fib.vs.athero.SMC <- log2(fib.vs.athero.SMC + 0.1)


## plot scatter
to.plot <- c("CDH1")
ggplot(data = fib.vs.athero.endo, aes(x = fib.vs.athero.endo[,1], y = fib.vs.athero.endo[,2])) +
  geom_point() + 
  geom_point(data = fib.vs.athero.endo[to.plot,1:2], aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]), col = "red", size = 4) +
  geom_text(aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]),label = as.character(to.plot), hjust=0, vjust=0, col = "red") +
  ggtitle("Fibrous vs Atheromatous endo") + 
  xlab("Fibrous") + 
  ylab("Atheromatous")

to.plot <- c('CDH1')
ggplot(data = fib.vs.athero.SMC, aes(x = fib.vs.athero.SMC[,1], y = fib.vs.athero.SMC[,2])) +
  geom_point() + 
  geom_point(data = fib.vs.athero.endo[to.plot,1:2], aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]), col = "red", size = 4) +
  geom_text(aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]),label = as.character(to.plot), hjust=0, vjust=0, col = "red") +
  ggtitle("Fibrous vs Atheromatous SMCs") + 
  xlab("Fibrous") + 
  ylab("Atheromatous")

#-----------------------------------------------------------------------------------------
# 2019-4-15 plot monocle

seuset <- readRDS("SMCs and ECs 18 patients.RDS")

TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "Endo and SMCs subset")
to.plot <- c("SNAI2", "CDH5", "MYH11", "MCAM", "ENG", "CD34")
for (i in to.plot) {
  FeaturePlot(seuset,
              features.plot = i)
}


DoHeatmap(seuset,
          genes.use = c("SNAI2", "CDH5", "MYH11", "MCAM", "ENG", "CD34"),
          title = "Endo and SMCs subset")

#-----------------------------------------------------------------------------------------
#senescence
# gH2AX, p21, laminB1, p53, pRB, pH3, CTGF, PAI1, p16. Ki67

to.plot <- c("H2AFX", "CDKN1A", "LMNB1", "TP53", "RB1", "PHC3", "CTGF", "SERPINE1", "CDKN2A", "MKI67")

pdf(paste(result.folder, "/", Sys.Date(), " Senescence Markers all clusters.pdf", sep = ""))
DoHeatmap(seuset,
          genes.use = to.plot,
          title = "Senescence markers all clusters")

TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "all clusters")

for (i in to.plot) {
  FeaturePlot(seuset,
              features.plot = i)
}

dev.off()
#-----------------------------------------------------------------------------------------
# intermediate cluster
seubset <- readRDS("SMCs and ECs 18 patients.RDS")

DEG <- read.table("DEG one vs all.txt")
DEG.pos <- DEG[DEG$p_val_adj < 0.05 & DEG$cluster == "Intermediate" & DEG$avg_logFC > 0,]
DEG.neg <- DEG[DEG$p_val_adj < 0.05 & DEG$cluster == "Intermediate" & DEG$avg_logFC < 0,]

plot.pos <- readClipboard()
plot.pos

pdf(paste(result.folder, "/", Sys.Date(), "pos DEG intermediate cells all communities.pdf", sep = ""))
lapply(plot.pos, function(x){
  VlnPlot(seuset,
          features.plot = x)
})
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), "pos DEG intermediate cells SMCS and EC subset.pdf", sep = ""))
lapply(plot.pos, function(x){
  VlnPlot(seubset,
          features.plot = x)
})
dev.off()

plot.neg <- readClipboard()
plot.neg 

pdf(paste(result.folder, "/", Sys.Date(), "neg DEG intermediate cells all communities.pdf", sep = ""))
lapply(plot.neg, function(x){
  VlnPlot(seuset,
          features.plot = x)
})
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), "neg DEG intermediate cells SMCS and EC subset.pdf", sep = ""))
lapply(plot.neg, function(x){
  VlnPlot(seubset,
          features.plot = x)
})
dev.off()
dev.off()


actins <- c("ACTA1", "ACTA2", "ACTB", "ACTBL2", "ACTC1", "ACTG1", "ACTG2")
lapply(actins, function(x){
  VlnPlot(seubset,
          features.plot = x)
})
dev.off()

#---------------------------------------------------------------------------------
# 9-5-2019

TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Patient", plot.title = "Patient")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "Ident")

f <- read.delim("marker genes (significant) 16 communities.txt")
View(f)

#-----------------------------------------------------------------------------------------

average.seuset <- AverageExpression(seuset, return.seurat = T)

endo <- rownames(seuset@meta.data[seuset@meta.data$res.1.7 == 10 | seuset@meta.data$res.1.7 ==13,])

pdf(paste(result.folder, "/", Sys.Date(), "THBS1 per patient endo cells only Lianne.pdf", sep = ""))
DoHeatmap(seuset, genes.use = "THBS1", group.by = "Patient", slim.col.label = T, cells.use = endo,
          group.label.rot = T, group.spacing = .5, title = "endo cells only")
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), "THBS1 per patient all cells Lianne.pdf", sep = ""))
DoHeatmap(seuset, genes.use = "THBS1", group.by = "Patient", slim.col.label = T, 
          group.label.rot = T, group.spacing = .5, title = "all cells")
dev.off()

#-----------------------------------------------------------------------------------------

seuset <- readRDS("SMCs and ECs 18 patients.RDS")

pdf(paste(result.folder, "/", Sys.Date(), "Plots for Michal Vln.pdf", sep = ""))
to.plot <- c("MYH11", "CD34", "ACTA2", "ACKR1", "COL6A1", "PECAM1", "CD59", "SPARC")
VlnPlot(seuset, features.plot = to.plot)

to.plot <- c("SNAI1", "SNAI2", "TWIST1", "ZEB1")
VlnPlot(seuset, features.plot = to.plot)

dev.off()

pdf(paste(result.folder, "/", Sys.Date(), "Plots for Michal feature plot1.pdf", sep = ""))
to.plot <- c("MYH11", "CD34", "ACTA2", "ACKR1", "COL6A1", "PECAM1", "CD59", "SPARC")
FeaturePlot(seuset, features.plot = to.plot)

to.plot <- c("SNAI1", "SNAI2", "TWIST1", "ZEB1")
FeaturePlot(seuset, features.plot = to.plot)

dev.off()

pdf(paste(result.folder, "/", Sys.Date(), "Plots for Michal feature plot2.pdf", sep = ""))
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "Ident")
to.plot <- c("MYH11", "CD34", "ACTA2", "ACKR1", "COL6A1", "PECAM1", "CD59", "SPARC", "SNAI1", "SNAI2", "TWIST1", "ZEB1")
for (i in to.plot) {
  FeaturePlot(seuset, features.plot = i)
}
dev.off()

singler <- CreateSinglerObject(counts = seuset@raw.data[,rownames(seuset@meta.data)],       
                               annot = seuset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "SMCs", 
                               min.genes = 0,                          # if starting from seurat, set min.genes to 0 
                               technology = "CEL-seq2", 
                               species = "Human", 
                               citation = "",
                               ref.list = list(), 
                               normalize.gene.length = FALSE, 
                               variable.genes = "de",
                               fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
                               do.signatures = TRUE,
                               clusters = seuset@ident,
                               do.main.types = TRUE, 
                               reduce.file.size = TRUE, 
                               numCores = 4)

pdf(paste(result.folder, "/", Sys.Date(), "Plots for Michal SingleR.pdf", sep = ""))
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 20)
dev.off()

#-----------------------------------------------------------------------
seubset <- readRDS("SMCs and ECs 18 patients.RDS")

to.plot <- c("CD34", "PECAM1", "TIE1", "VCAM1", "MYH11", "ACTA2")

pdf(paste(result.folder, "/", Sys.Date(), "Marie endoMT.pdf", sep = ""))
VlnPlot(seuset, features.plot = to.plot)
dev.off()

load(file = "all.seur.combined.E_clusters.rdata")
seuset <- all.seur.combined.E_clusters

singler <- CreateSinglerObject(counts = seuset@raw.data[,rownames(seuset@meta.data)],       
                               annot = seuset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "SMCs", 
                               min.genes = 0,                          # if starting from seurat, set min.genes to 0 
                               technology = "CEL-seq2", 
                               species = "Human", 
                               citation = "",
                               ref.list = list(), 
                               normalize.gene.length = FALSE, 
                               variable.genes = "de",
                               fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
                               do.signatures = TRUE,
                               clusters = seuset@ident,
                               do.main.types = TRUE, 
                               reduce.file.size = TRUE, 
                               numCores = 4)


VlnPlot(all.seur.combined, ident.include = "MYH11+ Smooth Muscle Cells", features.plot = to.plot)


to.plot <- c("CD34", "PECAM1", "CDH5", "BMP4", "ACKR1", "ENG", "VWF", "NOTCH1", "FGF18", "FN1", "MYH11", "ACTA2" )
pdf(paste(result.folder, "/", Sys.Date(), " Marie EAS plots.pdf", sep = ""))
VlnPlot(seuset, features.plot = to.plot)
dev.off()

#-----------------------------------------------------------------------
seuset <-  load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined

to.plot <- c("HDAC9", "IL6R", "PRKAR2B", "HBP1", "COG5", "GPR22", "DUS4L", "BCAP29", "SLC26A4", 
             "SLC26A4-AS1", "CBLL1", "SLC26A3", "DLD", "LAMB1" ,"GUCY1A1", "NOS3")

pdf(paste(result.folder, "/", Sys.Date(), " Voor Sander EAS.pdf", sep = ""))
lapply(to.plot, function(x){
  FeaturePlot(seuset, features.plot = x)
})
dev.off()

to.plot <- c("COL4A1", "COL4A2")

pdf(paste(result.folder, "/", Sys.Date(), " Voor Sander EAS2.pdf", sep = ""))
lapply(to.plot, function(x){
  FeaturePlot(seuset, features.plot = x)
})
dev.off()

#-----------------------------------------------------------------------
# 2019-6-3 Robin

seuset <-  load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined

average.seuset <- AverageExpression(seuset, return.seurat = T)
View(average.seuset@scale.data)

write.table(average.seuset@scale.data, "18 Patients average expression.txt")

#-----------------------------------------------------------------------
seuset <- readRDS("SMCs and ECs 18 patients.RDS")

pdf(paste(result.folder, "/", Sys.Date(), " EC and SMCs PECAM1, MYH11, ACTA2, CD34.pdf", sep = ""))
VlnPlot(seuset, features.plot = c("MYH11", "PECAM1", "ACTA2", "CD34"))
dev.off()

#-----------------------------------------------------------------------
# 2019-6-3 DEG endotheel 1 and endotheel 2

DEG <- FindMarkers(seuset, 
                   ident.1 = "CD34+ Endothelial Cells", 
                   ident.2 = "CD34+ Endothelial Cells II",
                   #test.use = "DESeq2",
                   logfc.threshold = 0.25,
                   min.pct = 0.1)

View(DEG)
DEG.pos <- subset(DEG, avg_logFC > 0)
DEG.neg <- subset(DEG, avg_logFC < 0)

write.table(DEG, "DEG endo1 en endo2 (Koen).txt")

#-----------------------------------------------------------------------------
# 2019-6-14

seuset <- readRDS("SMCs and ECs 18 patients.RDS")

pdf(paste(result.folder, "/", Sys.Date(), " EC markers for indermediate cluster.pdf", sep = ""))
VlnPlot(seuset, features.plot = c("CD34", "PECAM1", "VWF", "CDH5", "CLDN5", "NOS3"))
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), " SMC markers for indermediate cluster.pdf", sep = ""))
VlnPlot(seuset, features.plot = c("ACTA2", "MYH11"))
dev.off()

clip <- readClipboard()
clip2 <- readClipboard()
clip3 <- readClipboard()
clip4 <- readClipboard()

pdf(paste(result.folder, "/", Sys.Date(), " markers for intermediate cluster.pdf", sep = ""))
VlnPlot(seuset, features.plot = clip)
VlnPlot(seuset, features.plot = clip2)
# VlnPlot(seuset, features.plot = clip3)
# VlnPlot(seuset, features.plot = clip4)
for (i in clip){
FeaturePlot(seubset, features.plot = i)
}
for (i in clip2){
  FeaturePlot(seubset, features.plot = i)
}

# for (i in clip3){
#   FeaturePlot(seubset, features.plot = i)
# }
# for (i in clip4){
#   FeaturePlot(seubset, features.plot = i)
# }
"IL1RL1"

pdf(paste(result.folder, "/", Sys.Date(), " markers for intermediate cluster CD34 and IL1RL1.pdf", sep = ""))
FeaturePlot(seubset, features.plot = c("CD34"))
FeaturePlot(seubset, features.plot = c("IL1RL1"))
dev.off()

FeaturePlot(seubset, features.plot = "CD34")

#-----------------------------------------------------------------------------
## finding markers for paraffin staining of endoMT cluster

seuset <- readRDS("SMCs and ECs 18 patients.RDS")
seubset <- readRDS("seuset 18 patients 16 communities2.RDS")

pdf(paste(result.folder, "/", Sys.Date(), " Markers for paraffin staining endoMT.pdf", sep = ""))
DoHeatmap(seuset,
          genes.use = c("PECAM1", "ACTA2", "VWF", "MYH11", "IL1RL1"))

SMC <- c("NOTCH3", "EDNRB", "TAGLN")
VlnPlot(seuset, features.plot = SMC)

EC <- c("CDH5", "CDH1", "SELE", "SELP")
VlnPlot(seuset, features.plot = EC)

FeaturePlot(seubset, features.plot = c("ACTA2", "VWF", "MYH11", "IL1RL1"))
FeaturePlot(seubset, features.plot = "PECAM1")
dev.off()

#-----------------------------------------------------------------------------
# Voor Arjan

markers <- FindAllMarkers(all.seur.combined, logfc.threshold = 0.6, min.pct = 0.25)
write.table(markers, file = paste(result.folder, "/", Sys.Date(), " Marker genes all.seur.combined.txt", sep = ""))

#-----------------------------------------------------------------------------
## UMAP
load(file = "all.seur.combined.rdata")
scale.data <- as.matrix(all.seur.combined@scale.data)

umaps <- umap(t(scale.data), 
               n_neighbors = 57,
               scale = F,
               init = "pca",
               pca = 50)  #increases speed
plot(umaps) 

input.ggplot <- as.data.frame(unlist(umaps))
input.ggplot$cluster <- all.seur.combined@ident

pdf(paste(result.folder, "/", Sys.Date(), " UMAP all seur combined.pdf", sep = ""), width = 10)
ggplot() +
  geom_point(data = input.ggplot, aes(x = V1, y = V2, group = cluster, col = cluster), size = 0.4)
dev.off()

#-----------------------------------------------------------------------------

load(file = "all.seur.combined.rdata")

all.seur.combined <- AddMetaData(all.seur.combined, all.seur.combined@ident, col.name = "given.identity")

# CASP3


f <- as.matrix(all.seur.combined@raw.data)
f <- f[c("C3", "CASP3"),colnames(all.seur.combined@scale.data)]
dim(f)
g <- f[,f[1,] != 0 & f[2,] != 0]
dim(g)
double.pos <- colnames(g)
all.seur.combined@meta.data[double.pos,]

h <- f[,f[1,] != 0]
dim(h)
C3.pos <- colnames(h)
View(all.seur.combined@meta.data[C3.pos,])
View(table(all.seur.combined@meta.data[C3.pos,"given.identity"]))

j <- h[,h[1,] > 2]
View(table(all.seur.combined@meta.data[colnames(j),"given.identity"]))


pdf(paste(result.folder, "/", Sys.Date(), " C3 Nicholas Human data.pdf", sep = ""))
VlnPlot(all.seur.combined, features.plot = "C3", x.lab.rot = T)
TSNEPlot(all.seur.combined, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "ident")
FeaturePlot(all.seur.combined, features.plot = "C3")
FeaturePlot(all.seur.combined, features.plot = c("C3", "CASP3"), overlay = T)
dev.off()
