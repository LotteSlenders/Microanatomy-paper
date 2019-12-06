### 2019-4-18

# Project SCS 18 patients Experiment Guinea

# prepare seurat data object wit threshold nGene min 200

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Guinea")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.8
library(org.Hs.eg.db)   # v3.7.0
library(SingleR)        # v0.2.2

#-----------------------------------------------------------------------------------------
### prep meta data and files
# load meta data
meta.data <- read.table('Project SCS 18 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
# rename first column
names(meta.data)[1] <- "Patient"
# add ID column
meta.data$ID <- paste(meta.data$Patient, ".P",meta.data$Plate, sep = "")
# transform patient colum from integer to character
meta.data$Patient <- as.character(meta.data$Patient)

sum(duplicated(meta.data$ID))
View(meta.data)

# check files
all.files <- list.files(pattern="*.tsv")
all.files
all.files %in% meta.data$File
length(all.files) == dim(meta.data)[1]

#-----------------------------------------------------------------------------------------
### read in files and update gene names
## update symbols
#Set number of cells / plate.
n = 384
update.symbols <- function(my.counts.table = my.counts.table, n = n){
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping
  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(my.counts.table), columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- gene.info[!d,]
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- gene.info$SYMBOL
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  #Keep only highest expressed SYMBOL if multiple old aliases are now merged
  #Get maximum number of dups per gene
  num_dups <- max(table(my.counts.table$SYMBOL))
  
  #Loop until only one entry per gene left
  while(num_dups > 1){
    #First retrieve the list of duplicated symbols
    o <- order(my.counts.table$SYMBOL)
    my.counts.table <- my.counts.table[o,]
    d <- duplicated(my.counts.table$SYMBOL)
    
    #Then keep only the highest expressed one (this part will fail if we happen upon 2 rows with equal expression, but for now that hasn't happened yet so meh ;)
    h <- rowSums(my.counts.table[d,1:n]) > rowSums(my.counts.table[which(d==T)-1,1:n])
    h <- c(h, rowSums(my.counts.table[which(d==T)-1,1:n]) > rowSums(my.counts.table[d,1:n]))
    h <- h[which(!h)]
    my.counts.table <- my.counts.table[!row.names(my.counts.table) %in% names(h),]
    
    num_dups <- num_dups - 1 #One dup removed
  }
  
  #Overwrite old symbols in row names and clean up
  row.names(my.counts.table) <- my.counts.table$SYMBOL
  my.counts.table$SYMBOL <- NULL
  
  return(my.counts.table)
}

## read in and process
seurat.object.list <- list(c)
unwanted.genes <- c("^ERCC", "^MT-","^RPL", "^RPS","^UGDH-AS1$","^PGM2P2$", "^LOC100131257$", "^KCNQ1OT1$", 
                    "^MALAT1$", "^PGM5P2$", "^MAB21L3$","^EEF1A1$","^MALAT1$")

for (df in seq_along(all.files)) {
  current.df <- read.delim(all.files[df])
  
  # fix gene names and rownames
  current.df$GENEID <- sub("__.*", "" ,  current.df$GENEID)
  rownames(current.df) <- current.df$GENEID
  current.df <- current.df[,-1]
  # remove unwanted genes
  current.df <- current.df[grep(paste(unwanted.genes, collapse = "|"),rownames(current.df),invert=TRUE),]
  
  # match AE and plate number to file
  current.meta.data <- meta.data[meta.data$File == all.files[df],]
  colnames(current.df) <- paste(current.meta.data$ID, ".", 1:length(current.df), sep = "")
  
  # update gene names
  current.df <- update.symbols(current.df, n)
  
  ## create object
  current.object <- CreateSeuratObject(current.df, project = current.meta.data$ID)
  # add corresponding meta data
  current.object@meta.data <- cbind(current.object@meta.data, current.meta.data)
  # add to list
  seurat.object.list[[df]] <- current.object
  
  if (df == round(0.5*length(all.files))) {
    print("Halfway done")
  } 
  
  if (df == length(all.files)){
    print(":)")
  }
  
}

# somehow creating the objects prior occurs with a loss of cells for some, I do not know the cause
View(head(current.df))
head(current.object@meta.data)
tail(current.object@meta.data)
dim(current.object@meta.data)
dim(current.object@raw.data)
dim(current.object@data)

#-----------------------------------------------------------------------------------------
### QC
# merge suerat objects from list
seuset <- seurat.object.list[[1]]
for (i in 2:length(x = seurat.object.list)) {
  # only filter genes with too low expression during last merge
  if (i == length(seurat.object.list)){
    seuset <- MergeSeurat(object1 = seuset, object2 = seurat.object.list[[i]], do.normalize = F,
                          min.cells = 10)
  } else {
    seuset <- MergeSeurat(object1 = seuset, object2 = seurat.object.list[[i]], do.normalize = F)
  }
}
# remove object list to save memory
#rm(seurat.object.list)

dim(seuset@data)

## Filter cells with too low and too high gene counts
# plot nGene distribution
upper.bound <- mean(log2(seuset@meta.data$nGene)) + 2 * sd(log2(seuset@meta.data$nGene))
lower.bound <- mean(log2(seuset@meta.data$nGene)) - 2 * sd(log2(seuset@meta.data$nGene))

plot(density(log2(seuset@meta.data$nGene)))
abline(a = NULL, b = NULL, v = upper.bound, col = "red")
abline(a = NULL, b = NULL, v = lower.bound, col = "red")
abline(a = NULL, b = NULL, v = mean(log2(seuset@meta.data$nGene)), col = "blue")

2^upper.bound
2^lower.bound
mean(seuset@meta.data$nGene)

# vln plot
VlnPlot(object = seuset, features = c("nGene"))

# per plate
input.ggplot <- reshape2:: melt(seuset@meta.data[, c("nGene", "ID")])
graph.ggplot <- ggplot() +
  geom_density(input.ggplot, mapping = aes(group = ID, color = ID, x = value), stat = "density") +
  xlim(0, 750)
graph.ggplot
# interactive plot
#plotly::ggplotly(graph.ggplot) 

input.ggplot <- reshape2:: melt(seuset@meta.data[, c("nGene", "ID")])
input.ggplot$value <- log2(input.ggplot$value)
graph.ggplot <- ggplot() +
  geom_density(input.ggplot, mapping = aes(group = ID, color = ID, x = value), stat = "density") +
  xlim(0, 15)
graph.ggplot

# nUMI
plot(density(log2(seuset@meta.data$nUMI)))
abline(a = NULL, b = NULL, v = mean(log2(seuset@meta.data$nUMI)), col = "blue")

# check #cells
test.cells <- function(number){
  test <- subset(seuset@meta.data, nGene > number)
  print(table(test$Patient))
  print(dim(test))
}
test.cells(400)
test.cells(500)
test.cells(200)

# pick thresholds
seuset <- FilterCells(object = seuset,
                      subset.names = c("nGene"),
                      low.thresholds = 200,
                      high.thresholds = 10000
)
dim(seuset@data)

# Normalize data
seuset <- NormalizeData(object = seuset,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000
)

# scale data (regress out UMI)
seuset <- ScaleData(object = seuset, 
                    vars.to.regress = c("nUMI"))

# Find variable genes 
seuset <- FindVariableGenes(seuset, do.plot = F)
head(seuset@var.genes, n = 20)

#-----------------------------------------------------------------------------------------
### Clustering ICA
# # Run ICA
# random.seed <- set.seed(99)
# seuset <- RunICA(object = seuset,
#                  ics.compute = 50,
#                  genes.print = 1,
#                  seed.use = random.seed)

# seuset.ICA <- FindClusters(object = seuset, 
#                            reduction.type = "ica", 
#                            resolution = 1.7, 
#                            dims.use = 1:15, 
#                            force.recalc = T)
# 
# seuset.ICA <- RunTSNE(seuset.PCA,
#                       reduction.use = "ica",
#                       dims.use = 1:15
# )
# 
# TSNEPlot(seuset.ICA, do.label = T, pt.size = 1, label.size = 5, group.by = "ident")

#-----------------------------------------------------------------------------------------
### PCA reduction
random.seed <- set.seed(91)

# Run PCA
seuset <- RunPCA(object = seuset,
                 do.print = F)

PCElbowPlot(seuset, num.pc = 20)

# check PCA
PCAPlot(seuset,
        group.by = "Phenotype")
PCAPlot(seuset,
        group.by = "Patient")
PCAPlot(seuset,
        group.by = "Calcification")
PCAPlot(seuset,
        group.by = "Sex")

#-----------------------------------------------------------------------------------------
### clustering
seuset <- FindClusters(object = seuset, 
                       reduction.type = "pca", 
                       resolution = 1.7, 
                       dims.use = 1:10, 
                       force.recalc = T,
                       random.seed = 91 )

seuset <- RunTSNE(seuset,
                  reduction.use = "pca",
                  dims.use = 1:10,
                  seed.use = random.seed)

## check stuff
# how many cells per type
sort(table(seuset@meta.data$res.1.7))
# how many cells per plate
sort(table(seuset@meta.data$ID))
# how many cells per type per plate
table(seuset@meta.data$res.1.7, seuset@meta.data$ID)
# how many cells per type per patient
table(seuset@meta.data$res.1.7, seuset@meta.data$Patient)
# how many cells per patient
sort(table(seuset@meta.data$Patient))

barplot(prop.table(x = table(seuset@ident, seuset@meta.data$Patient)))
barplot(prop.table(x = table(seuset@ident, seuset@meta.data$ID)))

## Check tSNE 
pdf(paste(result.folder, "/", Sys.Date(), " Preliminary plots 18 communities.pdf", sep = ""))
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")


TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Calcification", plot.title = "Calcification")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Patient", plot.title = "AE number")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ID", plot.title = "AE and plate ID")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Phenotype", plot.title = "Plaque phenotype")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "SR.score", plot.title = "SR score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "CD68.score", plot.title = "CD68 score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "CD3.score", plot.title = "CD3 score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "CD34.score", plot.title = "CD34 score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "alpha.SMA.score", plot.title = "alpha-SMA score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Glyc.c.score", plot.title = "Glyc-c score")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Sex", plot.title = "Sex")

to.plot <- c("CD3E", "MYH11", "CD68", "CD14", "CD34", "KIT", "CD79A")
for (gene in to.plot) {
  FeaturePlot(seuset, features.plot = gene)
}

dev.off()

#-----------------------------------------------------------------------------------------
### infer celltypes
## single R
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
pdf(paste(result.folder, "/", Sys.Date(), " SingleR heatmap 18 communities.pdf", sep = ""))
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)
dev.off()

## Find marker genes
marker.genes <- FindAllMarkers(object = seuset,
                               logfc.threshold = 0.25,
                               min.pct = 0.1,
                               only.pos = T
)

#write.table(marker.genes, file = "marker genes 16 communities.txt", sep = "\t")

# save in misc slot
seuset@misc <- marker.genes

# subset only significant genes
marker.genes <- subset(marker.genes, p_val_adj <= 0.05)
write.table(marker.genes, file = "marker genes (significant) 16 communities.txt", sep = "\t")

# change ident

#-----------------------------------------------------------------------------------------
### if statisfacory, safe

saveRDS(seuset, file = "seuset 18 patients 16 communities2.RDS")


#-----------------------------------------------------------------------------------------




