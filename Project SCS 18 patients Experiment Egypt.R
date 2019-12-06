### 2019-4-10

# Project SCS 18 patients Experiment Egypt

# Lianne SMCs and ECs

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Egypt")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.8
library(SingleR)        # v0.2.2
library(pcaExplorer)    # v2.8.1

# object
seuset <- readRDS("seuset 18 patients 16 communities2.RDS")
fibrous <- readRDS("Fibrous patients only (18p).RDS")
atheromatous <- readRDS("Atheromatous patients only (18p).RDS")

#-----------------------------------------------------------------------------------------
# ### Subset the data
# # subset based on phenotype
# meta.data <- seuset@meta.data
# fibrous.cells <- rownames(meta.data[meta.data$Phenotype == "Fibrous",])
# 
# fibrous <- SubsetData(object = seuset,
#                       cells.use = fibrous.cells)
# 
# atheromatous.cells <- rownames(meta.data[meta.data$Phenotype == "Atheromatous" | meta.data$Phenotype == "Fibro-atheromatous",])
# 
# atheromatous <- SubsetData(object = seuset,
#                            cells.use = atheromatous.cells)
# 
# dim(fibrous@meta.data)
# dim(atheromatous@meta.data)
# 
# # downstream Functions 
# downstream.functions <- function(object) {
#   object <- FindVariableGenes(object, do.plot = F)
#   print(head(object@var.genes, n = 20))
#   
#   # Run PCA
#   object <- RunPCA(object = object,
#                    do.print = F)
#   
#   print(PCElbowPlot(object, num.pc = 20))
#   
#   return(object)
# }
# 
# 
# # Find clusters
# find.clusters <- function(object) { 
#   Random.Seed <- set.seed(53)
#   object <- FindClusters(object = object, 
#                          reduction.type = "pca", 
#                          resolution = 1, 
#                          dims.use = 1:10, 
#                          force.recalc = T,
#                          random.seed = 53 )
#   
#   
#   object <- RunTSNE(object,
#                     reduction.use = "pca",
#                     dims.use = 1:10,
#                     seed.use = Random.Seed)
#   
#   #return(object)
# }
# 
# 
# ### Fibrous
# ## PCA etc
# fibrous <- downstream.functions(fibrous)
# 
# ## find clusters
# fibrous <- find.clusters(fibrous)
# 
# singler <- CreateSinglerObject(counts = fibrous@raw.data[,rownames(fibrous@meta.data)],       
#                                annot = fibrous@ident,                                        # if is NULL it takes 10 clusters
#                                project.name = "Single Cell Plaques", 
#                                min.genes = 0,                          # if starting from seurat, set min.genes to 0 
#                                technology = "CEL-seq2", 
#                                species = "Human", 
#                                citation = "",
#                                ref.list = list(), 
#                                normalize.gene.length = FALSE, 
#                                variable.genes = "de",
#                                fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
#                                do.signatures = TRUE,
#                                clusters = fibrous@ident,
#                                do.main.types = TRUE, 
#                                reduce.file.size = TRUE, 
#                                numCores = 4
# )
# 
# SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50) + grid::grid.text("Fibrous", x = 0.8, y = 0.95)
# 
# dev.off()
# pdf(paste(result.folder, "/", Sys.Date(), " SingleR heatmap Fibrous.pdf", sep = ""))
# SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50) + grid::grid.text("Fibrous", x = 0.8, y = 0.95)
# dev.off()
# 
# to.plot <- c("MYH11", "CD34", "CD14", "CD68","CD3E","CD79A", "KIT")
# 
# pdf(paste(result.folder, "/", Sys.Date(), " tSNE plots Fibrous.pdf", sep = ""))
# TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "fibrous")
# for (i in to.plot) {
#   FeaturePlot(object = fibrous,
#               features.plot = i)
# }
# dev.off()
# 
# table(fibrous@meta.data$res.1)
# 
# ## find and save marker genes
# fibrous.markers <- FindAllMarkers(fibrous,
#                                   logfc.threshold = 0.25,
#                                   min.pct = 0.1)
# fibrous@misc <- fibrous.markers
# 
# 
# ## set ident
# new.ident <- c("T-cells", "T-cells", "Various", "T-cells", "T-cells", "Monocytes", "Monocytes", "Monocytes",
#                "Smooth-Muscle-cells", "Endothelial-cells", "B-cells", "Endothelial-cells", "Mast-cells")
# 
# for (i in 0:12) {
#   fibrous <- RenameIdent(object = fibrous, old.ident.name = i, new.ident.name = new.ident[i + 1])
# }
# 
# fibrous <- AddMetaData(fibrous, metadata = fibrous@ident, col.name = "given.identity")
# 
# pdf(paste(result.folder, "/", Sys.Date(), " Fibrous new names.pdf", sep = ""))
# TSNEPlot(fibrous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "fibrous new names")
# dev.off()
# 
# ## save
# saveRDS(fibrous, file = "Fibrous patients only (18p).RDS")
# 
# ### atheromatous
# ## PCA etc
# atheromatous <- downstream.functions(atheromatous)
# 
# ## find clusters
# atheromatous <- find.clusters(atheromatous)
# 
# singler <- CreateSinglerObject(counts = atheromatous@raw.data[,rownames(atheromatous@meta.data)],       
#                                annot = atheromatous@ident,                                        # if is NULL it takes 10 clusters
#                                project.name = "Single Cell Plaques", 
#                                min.genes = 0,                          # if starting from seurat, set min.genes to 0 
#                                technology = "CEL-seq2", 
#                                species = "Human", 
#                                citation = "",
#                                ref.list = list(), 
#                                normalize.gene.length = FALSE, 
#                                variable.genes = "de",
#                                fine.tune = FALSE, #Turn back on later (off for speedy debugging purposes)
#                                do.signatures = TRUE,
#                                clusters = atheromatous@ident,
#                                do.main.types = TRUE, 
#                                reduce.file.size = TRUE, 
#                                numCores = 4
# )
# 
# SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50) + grid::grid.text("Atheromatous", x = 0.8, y = 0.95)
# 
# dev.off()
# pdf(paste(result.folder, "/", Sys.Date(), " SingleR heatmap atheromatous.pdf", sep = ""))
# SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50) + grid::grid.text("Atheromatous", x = 0.8, y = 0.95)
# dev.off()
# 
# to.plot <- c("MYH11", "CD34", "CD14", "CD68","CD3E","CD79A", "KIT")
# 
# pdf(paste(result.folder, "/", Sys.Date(), " tSNE plots Atheromatous.pdf", sep = ""))
# TSNEPlot(atheromatous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "atheromatous")
# for (i in to.plot) {
#   FeaturePlot(object = atheromatous,
#               features.plot = i)
# }
# dev.off()
# 
# table(atheromatous@meta.data$res.1)
# 
# ## find and save marker genes
# atheromatous.markers <- FindAllMarkers(atheromatous,
#                                   logfc.threshold = 0.25,
#                                   min.pct = 0.1)
# atheromatous@misc <- atheromatous.markers
# 
# 
# ## set ident
# new.ident <- c("T-cells", "T-cells", "T-cells", "Smooth-Muscle-cells", "Various", "Monocytes", "Monocytes", 
#                "Monocytes", "Endothelial-cells", "B-cells", "Endothelial-cells", "Mast-cells")
# 
# for (i in 0:11) {
#   atheromatous <- RenameIdent(object = atheromatous, old.ident.name = i, new.ident.name = new.ident[i + 1])
# }
# 
# atheromatous <- AddMetaData(atheromatous, metadata = atheromatous@ident, col.name = "given.identity")
# 
# pdf(paste(result.folder, "/", Sys.Date(), " Atheromatous new names.pdf", sep = ""))
# TSNEPlot(atheromatous, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "atheromatous new names")
# dev.off()
# 
# # save
# saveRDS(atheromatous, file = "Atheromatous patients only (18p).RDS")


#-----------------------------------------------------------------------------------------
### fibrous vs atheromatous epithelium and smcs
# check
table(fibrous@meta.data$given.identity)
table(atheromatous@meta.data$given.identity)

## select cells of interest
# function
get.cells <- function(object, ident) {
  # transform object@raw.data into matrix
  current.data <- as.matrix(object@raw.data)
  # get all colums where the colnames correspond to the cells with chosen ident number
  current.matrix <- current.data[ ,rownames(object@meta.data[object@meta.data$given.identity == ident,])]
  return(current.matrix)
  # delete from environment
  rm(current.data, current.matrix)
}

# create dataset per celltype per phenotype
# since there are no repliates, we cannot use DESeq
fib.vs.athero.endo <- data.frame(endo.fibrous = rowSums(get.cells(fibrous, "Endothelial-cells")),
                            endo.atheromatous = rowSums(get.cells(atheromatous, "Endothelial-cells")))
                            

fib.vs.athero.SMC <- data.frame(smc.fibrous = rowSums(get.cells(fibrous, "Smooth-Muscle-cells")),
                                 smc.atheromatous = rowSums(get.cells(atheromatous, "Smooth-Muscle-cells")))

#write.table(fib.vs.athero.endo, "Fib vs Athero Endo not transformed.txt", sep = "\t", col.names = T, row.names = T)
#write.table(fib.vs.athero.SMC, "Fib vs Athero SMC not transformed.txt", sep = "\t", col.names = T, row.names = T)
  
# transform
fib.vs.athero.endo <- log2(fib.vs.athero.endo + 0.1)
fib.vs.athero.SMC <- log2(fib.vs.athero.SMC + 0.1)

## plot distribution
input.ggplot <- reshape2::melt(fib.vs.athero.endo, id.vars = NULL)
ggplot(data = input.ggplot , aes(x = variable, y = value)) +
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5, colour = "darkgrey") + 
  geom_boxplot(width=0.08, colour = "red") +
  xlab("phenotype") +
  ylab("log2 genecounts")

input.ggplot <- reshape2::melt(fib.vs.athero.SMC, id.vars = NULL)
ggplot(data = input.ggplot , aes(x = variable, y = value)) +
  geom_violin() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.5, colour = "darkgrey") + 
  geom_boxplot(width=0.08, colour = "red") +
  xlab("phenotype") +
  ylab("log2 genecounts")

## plot scatter
#to.plot <- c("SNAI2")
ggplot(data = fib.vs.athero.endo, aes(x = fib.vs.athero.endo[,1], y = fib.vs.athero.endo[,2])) +
  geom_point() + 
  #geom_point(data = fib.vs.athero.endo[to.plot,1:2], aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]), col = "red", size = 4) +
  #geom_text(aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]),label = as.character(to.plot), hjust=0, vjust=0, col = "red") +
  ggtitle("Fibrous vs Atheromatous endo") + 
  xlab("Fibrous") + 
  ylab("Atheromatous")

#to.plot <- c('SNAI2')
ggplot(data = fib.vs.athero.SMC, aes(x = fib.vs.athero.SMC[,1], y = fib.vs.athero.SMC[,2])) +
  geom_point() + 
  #geom_point(data = fib.vs.athero.endo[to.plot,1:2], aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]), col = "red", size = 4) +
  #geom_text(aes(x = fib.vs.athero.endo[to.plot,1], y = fib.vs.athero.endo[to.plot,2]),label = as.character(to.plot), hjust=0, vjust=0, col = "red") +
  ggtitle("Fibrous vs Atheromatous SMCs") + 
  xlab("Fibrous") + 
  ylab("Atheromatous")

#-----------------------------------------------------------------------------------------
### DESeq between fibrous patients and atheromatous patients
# prepare
input.meta.data <- read.table('Project SCS 18 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
names(input.meta.data)[1] <- "AE"
input.meta.data$ID <- paste(input.meta.data$AE, ".P",input.meta.data$Plate, sep = "")

# take only AE and phenotyope colum
patient.phenotype <- distinct(input.meta.data[,c(1,3,4)])
# if element in second colom equals x, replace with y, otherwise use original entry
patient.phenotype2 <- ifelse(patient.phenotype[,2] == "Fibro-atheromatous", "Atheromatous", patient.phenotype[,2])
# add to df
patient.phenotype$Phenotype2 <- patient.phenotype2
patient.phenotype$Sex.phenotype <- paste(patient.phenotype$Phenotype2, patient.phenotype$Sex, sep = ".")
patient.phenotype <- data.frame(patient.phenotype)
# remove value
rm(patient.phenotype2)  

# function to get the correct celltype per patient
get.cells <- function(object, AE.nr, identity) {
 
  # transform object@raw.data into matrix & get (raw) data
  current.matrix <- as.matrix(object@raw.data)
  #current.matrix <- as.matrix(object@data)
  # select cells that correspond to the correct ident and AE number
  current.cells <- subset(object@meta.data, object@meta.data$Patient == AE.nr)
  current.cells <- subset(current.cells, current.cells$given.identity == identity)
  current.cells <- rownames(current.cells)
  
  # print length current cells
  print(length(current.cells))
  
  # get the correct cells and sum 
  current.matrix <- rowSums(current.matrix[,current.cells, drop = F])
  
  return(current.matrix)
  
}

fib.vs.athero <- data.frame(matrix(nrow = dim(seuset@raw.data)[1], ncol = 0))
colnames.fib.vs.athero <- c()
#chosen.celltype <- "Endothelial-cells"
chosen.celltype <- "Smooth-Muscle-cells"
for (AE in patient.phenotype$AE) {
  # get phenotype
  if (patient.phenotype[patient.phenotype$AE == AE,2] == "Fibrous") {
    # get the cells
    correct.cells <- get.cells(fibrous, AE, chosen.celltype)
    name <- "fib"
  } else {
    correct.cells <- get.cells(atheromatous, AE, chosen.celltype)
    name <- "athero"
  }
  
  # bind to matrix
  fib.vs.athero <- cbind(fib.vs.athero, correct.cells)
  # save name
  colnames.fib.vs.athero <- c(colnames.fib.vs.athero, paste(AE,".",name, sep = ""))
  
  #remove from environment
  rm(correct.cells,AE, name)
}

# give colnames for clarity
names(fib.vs.athero) <- colnames.fib.vs.athero
head(fib.vs.athero)
# check
colSums(fib.vs.athero)
dim(fib.vs.athero)

## DESeq2 smcs and endo
# choose columns to proceed with
#chosen.patients <- c(2:7,9:14,16:18) # ECs
chosen.patients <- c(2:6, 10,11,14,17,18)  # SMCs
DESeq.fib.vs.athero <- fib.vs.athero[,chosen.patients]
head(DESeq.fib.vs.athero, n =1)
coldata <- patient.phenotype[chosen.patients,]

# round data
DESeq.fib.vs.athero <- round(DESeq.fib.vs.athero)

dds <- DESeqDataSetFromMatrix(countData =  DESeq.fib.vs.athero,
                              colData = coldata,
                              design = ~ Phenotype2)


dds <- DESeq(dds)                               
pre.results <- results(dds)
head(pre.results)
plotMA(pre.results)

results <- data.frame(pre.results)
results <- results[!is.na(results$pvalue),]

results <- results[order(results$pvalue),]
head(results)
sum(results$pvalue <= 0.05)

#write.table(results[results$pvalue <= 0.1,], "endo Fib vs Ath from raw.data.txt", sep = "\t", col.names = T, row.names = T)
#write.table(results[results$pvalue <= 0.1,], "SMCs Fib vs Ath from raw.data.txt", sep = "\t", col.names = T, row.names = T)


## PCA
transformed.dds <- vst(dds)
pcaplot(transformed.dds,
        intgroup = c("Phenotype2", "Sex"),
        ntop = 1000,
        text_labels = T,
        pcX = 1,
        pcY = 2)

pcaplot(transformed.dds,
        intgroup = c("Phenotype", "Sex"),
        ntop = 1000,
        text_labels = F,
        pcX = 1,
        pcY = 2)


#---------------------------------------------------------------------------------------------------
## DESeq sex
## DESeq2 smcs and endo
# choose columns to proceed with
#chosen.patients <- c(2:7,9:14,16:18) # ECs
chosen.patients <- c(2:6, 10,11,14,17,18)  # SMCs
DESeq.fib.vs.athero <- fib.vs.athero[,chosen.patients]
head(DESeq.fib.vs.athero, n =1)
coldata <- patient.phenotype[chosen.patients,]

# round data
DESeq.fib.vs.athero <- round(DESeq.fib.vs.athero)

dds <- DESeqDataSetFromMatrix(countData =  DESeq.fib.vs.athero,
                              colData = coldata,
                              design = ~ Sex.phenotype)


dds <- DESeq(dds)                               
pre.results <- results(object = dds,
                       contrast = c(1,3))
head(pre.results)
plotMA(pre.results)

results <- data.frame(pre.results)
results <- results[!is.na(results$pvalue),]

results <- results[order(results$pvalue),]
head(results)
sum(results$pvalue <= 0.05)
