### 2019-7-4

# Project SCS 18 patients Experiment Peru

# Monocle and slingshot on endoMT

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Peru")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(org.Hs.eg.db)   # v3.7.0
library(SingleR)        # v0.2.2
library(monocle)
library(slingshot)

# object
load(file = "all.seur.combined.rdata")

#-----------------------------------------------------------------------------------------
### get the cells
## Extract endo cell clusters
E_clusters <- c("CD34+ Endothelial Cells", "CD34+ Endothelial Cells II")
E_cells <- WhichCells(object = all.seur.combined, ident = E_clusters)

all.seur.combined.E_clusters <- SubsetData(object = all.seur.combined, ident.use = E_clusters)

# t-SNE and Clustering
all.seur.combined.E_clusters <- RunTSNE(all.seur.combined.E_clusters, reduction.use = "pca.aligned", dims.use = 1:15, 
                                        do.fast = T, perplexity = 17)
all.seur.combined.E_clusters <- FindClusters(all.seur.combined.E_clusters, reduction.type = "pca.aligned", 
                                             resolution = 0.8, dims.use = 1:15, force.recalc = T)

# set new ident
endo.ident <- c("ACKR1+ endothelial cells", "BMP4+ endothelial cells", "FGF18+ endothelial cells", "ACTA2+ endothelial cells")
for (i in 0:3) {
all.seur.combined.E_clusters <- RenameIdent(object = all.seur.combined.E_clusters, old.ident.name = i, new.ident.name = endo.ident[i + 1])
}

# new ident as meta data
all.seur.combined.E_clusters <- AddMetaData(all.seur.combined.E_clusters, all.seur.combined.E_clusters@ident, 
                                            col.name = "given.identity")
head(all.seur.combined.E_clusters@meta.data)

table(all.seur.combined.E_clusters@meta.data$given.identity, all.seur.combined.E_clusters@meta.data$Patient)

## extract smc clusters
SMC_clusters <- c("MYH11+ Smooth Muscle Cells")
SMC_cells <- WhichCells(object = all.seur.combined, ident = SMC_clusters)

all.seur.combined.SMC_clusters <- SubsetData(object = all.seur.combined, ident.use = SMC_clusters)


# t-SNE and Clustering
all.seur.combined.SMC_clusters <- RunTSNE(all.seur.combined.SMC_clusters, reduction.use = "pca.aligned", dims.use = 1:15, 
                                          do.fast = T)
all.seur.combined.SMC_clusters <- FindClusters(all.seur.combined.SMC_clusters, reduction.type = "pca.aligned", 
                                               resolution = 0.4, dims.use = 1:15, force.recalc = T)

# set new ident
smc.ident <- c("Synthetic smooth muscle cells", "Contractile smooth muscle cells")
for (i in 0:1) {
  all.seur.combined.SMC_clusters <- RenameIdent(object = all.seur.combined.SMC_clusters, old.ident.name = i, new.ident.name = smc.ident[i + 1])
}

# new ident as meta data
all.seur.combined.SMC_clusters <- AddMetaData(all.seur.combined.SMC_clusters, all.seur.combined.SMC_clusters@ident, 
                                            col.name = "given.identity")
head(all.seur.combined.SMC_clusters@meta.data)  

### merge objects
seuset <- MergeSeurat(all.seur.combined.E_clusters, all.seur.combined.SMC_clusters, do.scale = T)
head(seuset@meta.data)
tail(seuset@meta.data)

# ident
seuset <- SetIdent(seuset, ident.use = seuset@meta.data$given.identity)
head(seuset@ident)

## plot single R
# SingleR
singler <- CreateSinglerObject(counts = seuset@raw.data[,rownames(seuset@meta.data)],       
                               annot = seuset@ident,                                        # if is NULL it takes 10 clusters
                               project.name = "ECs and SMCs", 
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

SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 20)

pdf(paste(result.folder, "/", Sys.Date(), " singleR ECs and SMCs.pdf", sep = ""))
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 30)
dev.off()

table(seuset@meta.data$given.identity, seuset@meta.data$Patient)

#-----------------------------------------------------------------------------------------
### Monocle 
##create CellDataSet
# expression matrix
expression.matrix <- as.matrix(seuset@raw.data)
expression.matrix <- expression.matrix[,colnames(seuset@data)]
expression.matrix <- subset(expression.matrix, rowSums(expression.matrix) > 0) 

# phenodata
# remove unwanted variables
pd <- new("AnnotatedDataFrame", data = seuset@meta.data[colnames(expression.matrix), c(1,2,4,15)])
# featureData (one of the columns should be named "gene_short_name")
fd <- data.frame("gene_short_name" = row.names(expression.matrix),
                 row.names = row.names(expression.matrix))
fd <- new("AnnotatedDataFrame", data = fd)

# create dataset
dandy <- newCellDataSet(cellData = expression.matrix,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.1, #the minimum expression level that consistitutes true expression
                        expressionFamily = negbinomial()) # required 

# remove from environment
rm(pd,fd, expression.matrix)

### Estimate size factors and dispersions
dandy <- estimateSizeFactors(dandy)
# requires neg bin expression family. Is required for clustering
dandy <- estimateDispersions(dandy)

### filter
# detect genes
dandy <- detectGenes(dandy,
                     min_expr = 0.1)
head(fData(dandy))
max(fData(dandy)$num_cells_expressed)
min(fData(dandy)$num_cells_expressed) # can be low number because of subset
plot(density(fData(dandy)$num_cells_expressed))
dim(fData(dandy))

expressed_genes <- row.names(subset(fData(dandy),
                                    num_cells_expressed >= 5))

length(expressed_genes)
# plot num_cells_expressed
plot(density(fData(dandy)[expressed_genes,2]))

## filter cells out on number of genes expressed
# distribution of mRNA totals across the cells
# colsums total counts
pData(dandy)$Total_mRNAs <- Matrix::colSums(exprs(dandy))
# select only those with less than x total constructs
#dandy <- dandy[,pData(dandy)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(dandy)$Total_mRNAs)) +
                     2*sd(log10(pData(dandy)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(dandy)$Total_mRNAs)) -
                     2*sd(log10(pData(dandy)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(dandy), color = given.identity, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# remove cells with too low expression 
# argument not to do it, we can then directly project the states in seurat by adding them to meta.data
#dandy <- dandy[,pData(dandy)$Total_mRNAs > lower_bound & pData(dandy)$Total_mRNAs < 10000]
#dandy <- detectGenes(dandy, min_expr = 0.1)

# Once you've excluded cells that do not pass your quality control filters, 
#you should verify that the expression values stored in your CellDataSet follow a distribution that is roughly lognormal:
# log transform
L <- log(exprs(dandy[expressed_genes,])+ 0.1)
#L <- log(exprs(dandy[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- reshape2::melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "auto", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log") +
  ylab("Density")

### dimensional reduction
# reduce dimensions (DDRTree to avoid complications with pseudotime)
# overwrite dimensional reduction, but not clustering.
dandy <- reduceDimension(dandy, 
                         max_components = 2, 
                         num_dim = 5,
                         reduction_method = "DDRTree", 
                         residualModelFormulaStr = "~ nGene + Patient + nUMI + num_genes_expressed + num_genes_expressed + Total_mRNAs",
                         norm_method = "vstExprs",
                         verbose = F)

### pseudotime
## not cluster based
# step 1
# order cells otherwise next step does not work
dandy <- orderCells(dandy)
diff_test_res <- differentialGeneTest(dandy[expressed_genes,], verbose = T)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))                                    

# set them in the object
dandy <- setOrderingFilter(dandy, ordering_genes)
# plot
plot_ordering_genes(dandy)

# Step 2 reduce data dimensionality
# reduce space down to two dimensions
dandy <- reduceDimension(dandy, max_components = 2,
                         method = 'DDRTree')

# step 3 order cells along the trajectory
dandy <- orderCells(dandy)
#plot
#pdf(paste(result.folder, "/", Sys.Date(), " celltrajectory ECs and MCs.pdf", sep = ""))
#plot_cell_trajectory(dandy, color_by = "Cluster")
plot_cell_trajectory(dandy, color_by = "given.identity")
# after using orderCells, you can order bij cellState
plot_cell_trajectory(dandy, color_by = "State")
plot_cell_trajectory(dandy, color_by = "nGene")

dev.off()

# infer pseudotime
# identify the State which contact most of the cells from time zero, then pass this to orderCells
# no column with real 0 time point, unless there is a celltype (res) from which we can start.
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$res.0.8)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

dandy <- orderCells(dandy, root_state = GM_state(dandy))
# plot
plot_cell_trajectory(dandy, color_by = "Pseudotime")
# plot per state
plot_cell_trajectory(dandy, color_by = "State") +
  facet_wrap(~State, nrow = 1)

saveRDS(dandy, "monocle ECS and SMCs.RDS")


#-----------------------------------------------------------------------------------------
### Slingshot 
# http://www.bioconductor.org/packages/3.8/bioc/vignettes/slingshot/inst/doc/slingshot.html

# calculate PCA since it is not in the new object after combining
# var genes
seuset <- FindVariableGenes(seuset, do.plot = F)
head(seuset@var.genes, n = 20)

# Run PCA
seuset <- RunPCA(object = seuset,
                 do.print = F,
                 seed.use = 67,
                 pcs.compute = 15)

PCElbowPlot(seuset, num.pc = 20)

PCAPlot(seuset)

pca.coordinates <- seuset@dr$pca@cell.embeddings

# calculate tsne
seuset <- RunTSNE(seuset, reduction.use = "pca", dims.use = 1:15, 
                                          do.fast = T)
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "given identity")
TSNEPlot(all.seur.combined.E_clusters)
TSNEPlot(all.seur.combined.SMC_clusters)
#View(seuset@dr$tsne@cell.embeddings)

# # run UMAP
# reticulate::py_install(packages ='umap-learn')
# library(umap)
# seuset <- RunUMAP(seuset, dims.use = 1:15, reduction.use = "pca")

# gene filtering
expression.matrix <- as.matrix(seuset@raw.data)
expression.matrix <- expression.matrix[,colnames(seuset@data)]
expression.matrix <- subset(expression.matrix, rowSums(expression.matrix) > 3) 

View(seuset)

# both don't work
ammo <- slingshot(t(expression.matrix), reduceDim = "PCA", clusterLabels = seuset@meta.data$given.identity)
ammo <- slingshot(t(expression.matrix), reduceDim = "PCA")
summary(ammo$slingPseudotime_1)

#-----------------------------------------------------------------------------------------
### Diffusionmap
# https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html#diffusion-map-pseudotime
library(destiny)

# prepare count 
expression.matrix <- as.matrix(seuset@raw.data)
expression.matrix <- expression.matrix[,colnames(seuset@data)]
expression.matrix <- subset(expression.matrix, rowSums(expression.matrix) > 3)
expression.matrix <- log2(expression.matrix + 0.1)
head(expression.matrix)

# find optimal sigma
sigmas <- find_sigmas(t(expression.matrix), verbose = FALSE)
dm <- DiffusionMap(t(expression.matrix), sigma = optimal_sigma(sigmas))

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = seuset@meta.data$given.identity)

ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

pseudotime.diffusionmap <- rank(eigenvectors(dm)[,1])
head(pseudotime.diffusionmap)
