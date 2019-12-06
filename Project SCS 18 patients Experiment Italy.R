### 2019-3-28

# Project SCS 18 patients Experiment Italy

# resolution check for SMCs and ECs subset without that one patient.

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Italy")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(SingleR)        # v0.2.2
library(monocle)        # v2.10.1

# object
seuset <- readRDS("SMCs and ECs 18 patients.RDS")

#-----------------------------------------------------------------------------------------
### get SMCs
# original data
seuset <- readRDS("seuset 18 patients 16 communities2.RDS")

# subset
seubset <- SubsetData(seuset,
                      ident.use = c(10, 13, 8))

# remove that one patient
table(seubset@meta.data$Patient, seubset@meta.data$res.1.7)

seubset <- SubsetData(seubset,
                      cells.use = rownames(seubset@meta.data[seubset@meta.data$Patient != 4443,]))

# find variable genes
seubset <- FindVariableGenes(seubset, do.plot = F)
# Run PCA
seubset <- RunPCA(seubset, do.print = F)
# Find clusters
seubset <- FindClusters(object = seubset,
                        reduction.type = "pca",
                        resolution = 1,
                        dims.use = 1:10,
                        force.recalc = T,
                        random.seed = 91 )
# run tSNE
seubset <- RunTSNE(seubset, dims.use = 1:10)



# SingleR
singler <- CreateSinglerObject(counts = seubset@raw.data[,rownames(seubset@meta.data)],       
                               annot = seubset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "SMCs and ECs", 
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

# plot
pdf(paste(result.folder, "/", Sys.Date()," ", "SMCs and ECs without patient 4443.pdf", sep = ""))

TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")
TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "Patient", plot.title = "Patient")
TSNEPlot(seubset, do.label = F, pt.size = 1, label.size = 5, group.by = "Phenotype", plot.title = "Phenotype")

to.plot <- c("CD3E", "MYH11", "CD68", "CD14", "CD34")
for (gene in to.plot) {
  FeaturePlot(seubset, features.plot = gene)
}

SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)

dev.off()
