### 2019-3-28

# Project SCS 18 patients Experiment Japan

# subset t cells for Marie

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Japan")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(SingleR)        # v0.2.2

# object
seuset <- readRDS("seuset 18 patients 16 communities2.RDS")

#-----------------------------------------------------------------------------------------
### get T cells
# set seed
set.seed(78)

# subset
seubset <- SubsetData(seuset,
                      ident.use = c(0,1,3,5,9))

# find variable genes
seubset <- FindVariableGenes(seubset, do.plot = F)

# Run PCA
seubset <- RunPCA(seubset, do.print = F)

# Find clusters
seubset <- FindClusters(object = seubset,
                        reduction.type = "pca",
                        resolution = 0.6,
                        dims.use = 1:10,
                        force.recalc = T,
                        random.seed = 91 )
# run tSNE
seubset <- RunTSNE(seubset, dims.use = 1:10)



# SingleR
singler <- CreateSinglerObject(counts = seubset@raw.data[,rownames(seubset@meta.data)],       
                               annot = seubset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "T cells", 
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
pdf(paste(result.folder, "/", Sys.Date(), " Tcells SCS 18 patients.pdf", sep = ""))

TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")
to.plot <- c("CD3E", "MYH11", "CD68", "CD14", "CD34")
for (gene in to.plot) {
  FeaturePlot(seubset, features.plot = gene)
}
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)

dev.off()

#-----------------------------------------------------------------------------------------
### add names for analysis
cluster.number <- 8

## set ident
new.ident <- c("Tcell00", "Tcell01", "Tcell02", "Tcell03", "Tcell04", "Tcell05", "Tcell06", "Tcell07")

for (i in 0:7) {
  seubset <- RenameIdent(object = seubset, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

# add meta data
seubset <- AddMetaData(seubset, metadata = seubset@ident, col.name = "given.identity")

#-----------------------------------------------------------------------------------------
### DEG
## one vs all
all.markers <- FindAllMarkers(seubset,
                              logfc.threshold = 0.6,
                              min.pct = 0.1 )

all.markers <- subset(all.markers, p_val_adj <= 0.05)

# save
write.table(all.markers, "one vs all T cells 18 patients.txt", sep = " ")

## one vs one
# prepare all possible combinations
choices <- choose(cluster.number, 2)        # identify the number of possible pairs
combinations <- combn(cluster.number, 2)    # create those pairs

differential.expression.clusters <- list()

# pairwise comparison for differential expression
for (choice in 1:choices) {
  # select the cluster numbers to run
  cluster1 <- combinations[,choice][1] # the first index of the combination
  cluster2 <- combinations[,choice][2] # the second index of the combination
  
  # input function must be character (if stated so)
  cluster1 <- new.ident[cluster1]
  cluster2 <- new.ident[cluster2]
  
  # run formula
  pairwise.comparison <- FindMarkers(object = seubset,
                                     ident.1 = cluster1 ,
                                     ident.2 = cluster2,
                                     logfc.threshold = 0.6, #default 0.25
                                     min.pct = 0.1,
                                     test.use = "wilcox") # default wilcox
  
  # add gene names for later
  pairwise.comparison$gene <- rownames(pairwise.comparison)
  # add clusters
  pairwise.comparison$cluster1 <- cluster1
  pairwise.comparison$cluster2 <- cluster2
  # select only values with p adjusted of certain value
  pairwise.comparison <- pairwise.comparison[(pairwise.comparison$p_val_adj <= .05), ]
  # add to list
  differential.expression.clusters[[choice]] <- pairwise.comparison
  
}

# save lists
lapply(differential.expression.clusters, function(x){
  # get the clusters
  cluster1 <- x$cluster1[1]
  cluster2 <- x$cluster2[2]
  
  # order by foldchange
  x <- x[order(x$avg_logFC, decreasing = T),]
  
  write.table(x, paste(cluster1, "vs", cluster2, ".txt", sep = " "))
})


to.plot <- c("FOXP3", "TBX21", "GATA3", "RORC", "RORA", "IFNG", 
             "ZNF683", "EOMES", "CTLA4", "LAG3", "HAVCR2", "PDCD1", "TIGIT", "ITGAE")

pdf(paste(result.folder, "/", Sys.Date(), " Markers Marie.pdf", sep = ""))
lapply(to.plot, function(x){
  VlnPlot(seubset, features.plot = x)
})

dev.off()


