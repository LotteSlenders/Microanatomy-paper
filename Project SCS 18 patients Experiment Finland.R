### 2019 - 4 - 24

# Project SCS 18 patients Experiment Finland

# SMCS and ECs specific genes: DEG

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Finland")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(SingleR)        # v0.2.2
library(monocle)        # v2.10.1

# object
seuset <- readRDS("SMCs and ECs 18 patients.RDS")

#-----------------------------------------------------------------------------------------
### prepare
cluster.number <- 5

## set ident
new.ident <- c("Endothelium2", "Smooth.muscle1", "Endothelium1", "Smooth.muscle2", "Intermediate")

for (i in 0:4) {
  seuset <- RenameIdent(object = seuset, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

# add meta data
seuset <- AddMetaData(seuset, metadata = seuset@ident, col.name = "given.identity")

# plot
#pdf(paste(result.folder, "/", Sys.Date(), " SMCs and ECs tSNE and singleR.pdf", sep = ""))
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "new ident")

singler <- CreateSinglerObject(counts = seuset@raw.data[,rownames(seuset@meta.data)],       
                               annot = seuset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "Single Cell Plaques", 
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
                               numCores = 4
)

SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)
dev.off()

#-----------------------------------------------------------------------------------------
### DEG
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
  
  # input function must be character
  cluster1 <- new.ident[cluster1]
  cluster2 <- new.ident[cluster2]
  
  # run formula
  pairwise.comparison <- FindMarkers(object = seuset,
                                     ident.1 = cluster1 ,
                                     ident.2 = cluster2,
                                     logfc.threshold = 0.6, #default 0.25
                                     min.pct = 0.2,
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

## one vs all
# differential expression findAllMarkers
differential.expression.all <- FindAllMarkers(object = seuset,
                                              logfc.threshold = 0.6,
                                              min.pct = 0.2,
                                              only.pos = FALSE,
                                              test.use = "wilcox")

dim(differential.expression.all)
differential.expression.all <- differential.expression.all[(differential.expression.all$p_val_adj <= 0.05),]
dim(differential.expression.all)

write.table(differential.expression.all, "DEG one vs all.txt", sep = " ")

# subset results for intermediate population
intermediate.DEG <- subset(differential.expression.all, differential.expression.all$cluster == "Intermediate")
intermediate.DEG.positive <- intermediate.DEG[intermediate.DEG$avg_logFC > 1, "gene"]
intermediate.DEG.negative <- intermediate.DEG[intermediate.DEG$avg_logFC < 1, "gene"]

PCAPlot(seuset,
        group.by = "ident")

DoHeatmap(seuset,
          genes.use = c("H2AFX", "CDKN1A", "LMNB1", "TP53", "RB1", "PHC3", "CTGF", "SERPINE1", "CDKN2A", "MKI67"),
          slim.col.label = T,
          group.by = "Patient",
          # cells.use = NULL,
          title = "titel",
          group.label.rot = T)
dev.off()          

