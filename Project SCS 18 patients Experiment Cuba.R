### 2019-3-28

# Project SCS 18 patients Experiment Cuba

# Monocle

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Cuba")

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
#pdf(paste(result.folder, "/", Sys.Date(), "Monocle SMCs and ECs 18 patients.pdf", sep = ""))

TSNEPlot(seubset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "identity")
to.plot <- c("CD3E", "MYH11", "CD68", "CD14", "CD34")
for (gene in to.plot) {
  FeaturePlot(seubset, features.plot = gene)
}
SingleR.DrawHeatmap(SingleR = singler$singler[[1]]$SingleR.clusters.main, clusters = levels(singler$meta.data$orig.ident), top.n = 50)

#dev.off()


# save subset
saveRDS(seubset, "SMCs and ECs 18 patients.RDS")

# remove from enviornment
rm(seuset)
# change
seuset <- seubset
# remove
rm(seubset)

#-----------------------------------------------------------------------------------------
### create CellDataSet
# object
seuset <- readRDS("SMCs and ECs 18 patients.RDS")

# expression matrix
expression.matrix <- as.matrix(seuset@raw.data)
expression.matrix <- expression.matrix[,colnames(seuset@data)]
expression.matrix <- subset(expression.matrix, rowSums(expression.matrix) > 0) 

# phenodata
# remove unwanted variables
pd <- new("AnnotatedDataFrame", data = seuset@meta.data[colnames(expression.matrix), c(1,2,4,6:13,17,19,20)])
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
min(fData(dandy)$num_cells_expressed)
plot(density(fData(dandy)$num_cells_expressed))
dim(fData(dandy))

expressed_genes <- row.names(subset(fData(dandy),
                                    num_cells_expressed >= 10))
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

qplot(Total_mRNAs, data = pData(dandy), color = Phenotype, geom =
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


### Classifying and Counting Cells
## Clustering cells without marker genes
# which genes to use?
disp_table <- dispersionTable(dandy)
# which genes to use for unsupervised clustering
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
dandy <- setOrderingFilter(dandy, unsup_clustering_genes$gene_id)
plot_ordering_genes(dandy)

### cluster
# clustering with DDRTree is not recommended. So for actual clustering use TSNE as reduction method. 
# determine PC number
plot_pc_variance_explained(dandy, return_all = F,
                           norm_method = "vstExprs",
                           max_components = 20) 

# for proper clusters
dandy <- reduceDimension(dandy, 
                         max_components = 2, 
                         num_dim = 6,
                         reduction_method = "tSNE", 
                         # residualModelFormulaStr = "~ nGene + nUMI ",
                         norm_method = "vstExprs",
                         verbose = F)

# cluster the cells
dandy <- clusterCells(dandy)
# plot results
plot_cell_clusters(cds = dandy,
                   x = 1, y = 2, 
                   color_by  = "res.1"
)

plot_cell_clusters(cds = dandy,
                   x = 1, y = 2, 
                   #markers = c(),
                   color_by  = "Cluster"
)


# reduce dimensions (DDRTree to avoid complications with pseudotime)
# overwrite dimensional reduction, but not clustering.
dandy <- reduceDimension(dandy, 
                         max_components = 2, 
                         num_dim = 6,
                         reduction_method = "DDRTree", 
                         # residualModelFormulaStr = "~ nGene + nUMI ",
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
#pdf(paste(result.folder, "/", Sys.Date(), " celltrajectory SMCs, ECs and MC.pdf", sep = ""))
plot_cell_trajectory(dandy, color_by = "Cluster")
plot_cell_trajectory(dandy, color_by = "res.1")
# after using orderCells, you can order bij cellState
plot_cell_trajectory(dandy, color_by = "State")
plot_cell_trajectory(dandy, color_by = "nGene")

dev.off()

# infer pseudotime
# identify the State which contact most of the cells from time zero, then pass this to orderCells
# no column with real 0 time point, unless there is a celltype (res) from which we can start.
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$res.1)[,"2"]
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

write.csv(diff_test_res, file = "monocle diff expr results.txt")
saveRDS(dandy, "monocle SMCs, ECs and MC subset 9 patients.RDS")

#-------------------------------------------------------------------------------------------------
### from monocle back to seurat
## add pData columns back to seurat 

seuset <- AddMetaData(object = seuset,
                      metadata = data.frame(pData(dandy)[,5:15]))
head(seuset@meta.data)

saveRDS(seuset, "Colon plus monocle.RDS")

pdf(paste(result.folder, "/", Sys.Date(), " monocle and seurat.pdf", sep = ""))
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "State", plot.title = "Monocle state")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "Cluster", plot.title = "Monocle clusters")
TSNEPlot(seuset, do.label = T, pt.size = 1, label.size = 5, group.by = "ident", plot.title = "ident")

dev.off()

