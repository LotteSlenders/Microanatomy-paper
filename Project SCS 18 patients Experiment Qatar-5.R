## 2019-10-21

# Project SCS 18 patients Experiment Qatar-5

# receptor ligand expression in single cell data
# DOI: 10.1038/ncomms8866
# subcluster analysis

# final version

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Qatar")

# create folder
result.folder <- paste(Sys.Date(), "results Quatar5", sep = " ")
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
#load(file = "all.seur.combined.rdata")
seuset <- readRDS("all.seur.split.R")
#seuset <- all.seur.combined


#-----------------------------------------------------------------------------------------
# ### change names for clarity
# # get the old names to change to new names
# celltypes <- c("CD68+ Myeloid cells II", "CD68+ Myeloid cells I", "CD68+ Myeloid cells III", "CD34+ Endothelial cells II",
#                "CD3+ T-Cells III", "CD68+ Myeloid cells IV", "CD3+ T-Cells II", "Mixed cells", "CD3+ T-Cells I" ,
#                "ACTA2+ Smooth Muscle Cells", "KIT+ Mast-Cells", "CD34+ Endothelial cells I", "CD3+ T-Cells VI", "CD79A+ B-cells")
# 
# for (i in as.character(unique(seuset@ident))) {
#   count <- which(as.character(unique(seuset@ident)) == i)
#   seuset <- RenameIdent(object = seuset, old.ident.name = i, new.ident.name = celltypes[count])
# }
# 
# unique(seuset@ident)
# 
# ### break up clusters
# # change ident function
# change.ident <- function(Object, new.ident.string){
#   for(i in unique(Object@ident)){
#     Object <- RenameIdent(object = Object, old.ident.name = i, new.ident.name = new.ident.string[as.numeric(i) + 1])
#   }
#   df <- data.frame(Object@ident)
#   return(df)
# }
# 
# CD4 <- change.ident(all.seur.combined.T_split_clusters.CD4, c("CD4.0", "CD4.1", "CD4.2", "CD4.3", "CD4.4"))
# CD8 <- change.ident(all.seur.combined.T_split_clusters.CD8, c("CD8.0", "CD8.1", "CD8.2", "CD8.3"))
# My <- change.ident(all.seur.combined.M_clusters, c("My.0", "My.1", "My.2", "My.3", "My.4"))
# SMC <- change.ident(all.seur.combined.SMC_clusters, c("S.0", "S.1"))
# EC <- change.ident(all.seur.combined.E_clusters, c("E.0", "E.1", "E.2","E.3"))
# 
# BC <- data.frame(Object.ident = subset(all.seur.combined@ident, all.seur.combined@ident == "CD79A+ B-Cells"))
# Mast <- data.frame(Object.ident = subset(all.seur.combined@ident, all.seur.combined@ident == "KIT+ Mast Cells"))
# Mixed <- data.frame(Object.ident = subset(all.seur.combined@ident, all.seur.combined@ident == "Mixed Cells"))
# 
# overview.cells <- rbind(CD4, CD8, My, SMC, EC, BC, Mast, Mixed)
# 
# ### make new seurat
# # we miss some of the t cells now so careful
# seuset <- SubsetData(all.seur.combined, cells.use = rownames(overview.cells))
# ## change ident
# # save old ident
# seuset <- AddMetaData(seuset, metadata = seuset@ident, col.name = "original.ident")
# # order overview.cells to match order seuset@ident
# overview.cells <- overview.cells[rownames(seuset@meta.data),,drop = F]
# # set new identity
# seuset <- SetIdent(seuset, cells.use = rownames(overview.cells), ident.use = overview.cells$Object.ident)
# # save new ident
# seuset <- AddMetaData(seuset, metadata = seuset@ident, col.name = "new.ident")
# 
# ## save long new idents 
# celltypes <- as.character(unique(overview.cells$Object.ident))
# long.ident <- c("IL1B+ inflammatory macrophages", "TNF+ inflammatory macrophages", "ABCG1+ foamy macrophages", 
#                 "CD1C+ dentritic cells", "CD3+CD4+ T-cells",  "BMP4+ endothelial cells", "FGF18+ endothelial cells",  
#                 "PRF1+ Cytotoxic CD4+ T cells", "IL7R+ Naive CD4+ T cells", "GZMB+ Cytotoxic CD8+ T cells",
#                 "CD8+ Effector/Memory T cells", "FOXP3+ Regulatory CD4+ T cells", "LEF1+ Central memory CD4+ T cells",
#                 "GZMK+ Cytotoxic CD4+ T cells", "IL7R+ Naive/Central MemoryCD8+ T cells", "Mixed cells", 
#                 "GZMK+ Effector/Memory CD8+ T cells", "Contractile smooth muscle cells", "Synthetic smooth muscle cells",
#                 "KIT+ Mast Cells", "ACKR1+ endothelial cells", "ACTA2+ endothelial cells", "CD79A+ B-Cells"
# )
# #g <- overview.cells
# overview.cells <- g
# 
# # 
# # 
# # for(i in 1:(length(celltypes))){
# #   j <- celltypes[i]
# #   print(j)
# #   overview.cells[j,"long.ident"] <- long.ident[i]
# # }
# 
# ### save object
# saveRDS(seuset, "all.seur.split.R")

#-----------------------------------------------------------------------------------------
### get the receptor ligand pairs
data <- read_excel("ncomms8866-s3 receptor ligand pairs.xlsx", sheet = "All.Pairs")
#View(data)
pairs <- data.frame(receptor = data$Receptor.ApprovedSymbol,
                    ligand = data$Ligand.ApprovedSymbol,
                    row.names = data$Pair.Name)
#View(pairs)

# check
dim(pairs)
sum(pairs$receptor %in% rownames(seuset@data))
sum(pairs$ligand %in% rownames(seuset@data))

# report
report <- data.frame(receptors_expressed = sum(pairs$receptor %in% rownames(seuset@data)),
                     ligands_epressed = sum(pairs$ligand %in% rownames(seuset@data)))

### remove items that are not present
pairs <- pairs[pairs$receptor %in% rownames(seuset@data) & pairs$ligand %in% rownames(seuset@data),]
dim(pairs)

# report
report$both_expressed <- dim(pairs)[1]
report$unique_receptors <- length(unique(pairs$receptor))
report$unique_ligands <- length(unique(pairs$ligand))
report

# write to file

#write.xlsx(report, file = paste(Sys.Date(), "Report_Q3.xlsx", sep = " "), sheetName = "all_pairs")

#-----------------------------------------------------------------------------------------
### create pseudobulk
# get celltypes to pseudobulk
#celltypes <- as.character(celltypes)
celltypes <- as.character(unique(seuset@ident))
celltypes

### get cells corresponding to celltype
# set minimum gene expression requirement
minimum.gene.expression <- 5
# get counts
gene.expression.matrix <- as.matrix(seuset@data)

pseudobulk.genes <- lapply(celltypes, function(x){
  # get cells corresponding to celltype
  cells <- WhichCells(seuset, ident = x)
  # select the cells from seuset@data
  my.matrix <- gene.expression.matrix[,cells]
  
  # set expression threshold in at least x% percent of cells
  threshold <- round(minimum.gene.expression/100 * length(cells))
  #print(threshold)
  
  # select genes that meet minimum gene expression requirement
  #my.matrix <- my.matrix[apply(my.matrix !=0, 1,sum) >= threshold,] # or
  my.matrix <- my.matrix[rowSums(my.matrix !=0) >= threshold,]
  
  # return selected genes for pseudobulk
  return(rownames(my.matrix))
  #return(my.matrix)
})

names(pseudobulk.genes) <- celltypes
# View(pseudobulk.genes)

#-----------------------------------------------------------------------------------------
### match per celltype the expression of receptor or ligand
# check the absolute numbers

expression.per.celltype <- lapply(celltypes, function(x){
  # get the genes that are expressed per celltype
  current.genes <- unlist(pseudobulk.genes[x])
  #print(dim(current.genes))
  
  # get the numbers
  df <- data.frame(receptor = sum(pairs$receptor %in% current.genes),
                   ligand = sum(pairs$ligand %in% current.genes),
                   both = dim(pairs[pairs$receptor %in% current.genes & pairs$ligand %in% current.genes,])[1],
                   unique_receptor = length(unique(pairs[pairs$receptor %in% current.genes,1])),
                   unique_ligand = length(unique(pairs[pairs$ligand %in% current.genes,2])),
                   total_genes = length(current.genes)
  )
})

names(expression.per.celltype) <- celltypes
#View(expression.per.celltype)


# save results
numbers.report <- do.call(rbind.data.frame, expression.per.celltype)
#write.xlsx(numbers.report, file = paste(Sys.Date(), "Report_Q3.xlsx", sep = " "), sheetName = "Per_celltype", append = T)

#-----------------------------------------------------------------------------------------
### which receptors and ligands are present

interaction.matrix2 <- lapply(celltypes, function(x){
  # get current celltype
  current.celltype <- unlist(pseudobulk.genes[x])
  
  # check if gene is present from pair
  current.interactions <- sapply(rownames(pairs), function(i){
    # select the current receptor and current ligand and check presence
    receptor.expression <- pairs[i,"receptor"] %in% current.celltype
    ligand.expression <- pairs[i, "ligand"] %in% current.celltype
    
    # bind
    table <- cbind(receptor.expression, ligand.expression)
    names(table) <- c("receptor", "ligand")
    
    return(table)
    
  })
  
  return(current.interactions)
  
})

names(interaction.matrix2) <- celltypes
#View(interaction.matrix2)
#View(interaction.matrix2[[1]])

#-----------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------
### match per celltype the expression of receptor or ligand
# check the absolute numbers

expression.per.celltype <- lapply(celltypes, function(x){
  # get the genes that are expressed per celltype
  current.genes <- unlist(pseudobulk.genes[x])
  #print(dim(current.genes))
  
  # get the numbers
  df <- data.frame(receptor = sum(pairs$receptor %in% current.genes),
                   ligand = sum(pairs$ligand %in% current.genes),
                   both = dim(pairs[pairs$receptor %in% current.genes & pairs$ligand %in% current.genes,])[1],
                   unique_receptor = length(unique(pairs[pairs$receptor %in% current.genes,1])),
                   unique_ligand = length(unique(pairs[pairs$ligand %in% current.genes,2])),
                   total_genes = length(current.genes)
  )
})

names(expression.per.celltype) <- celltypes
#View(expression.per.celltype)


# save results
numbers.report <- do.call(rbind.data.frame, expression.per.celltype)
#write.xlsx(numbers.report, file = paste(Sys.Date(), "Report_Q3.xlsx", sep = " "), sheetName = "Per_celltype", append = T)

#-----------------------------------------------------------------------------------------
### which receptors and ligands are present

interaction.matrix2 <- lapply(celltypes, function(x){
  # get current celltype
  current.celltype <- unlist(pseudobulk.genes[x])
  
  # check if gene is present from pair
  current.interactions <- sapply(rownames(pairs), function(i){
    # select the current receptor and current ligand and check presence
    receptor.expression <- pairs[i,"receptor"] %in% current.celltype
    ligand.expression <- pairs[i, "ligand"] %in% current.celltype
    
    # bind
    table <- cbind(receptor.expression, ligand.expression)
    names(table) <- c("receptor", "ligand")
    
    return(table)
    
  })
  
  return(current.interactions)
  
})

names(interaction.matrix2) <- celltypes
#View(interaction.matrix2)
#View(interaction.matrix2[[1]])

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
### match per celltype the expression of receptor or ligand
# check the absolute numbers

expression.per.celltype <- lapply(celltypes, function(x){
  # get the genes that are expressed per celltype
  current.genes <- unlist(pseudobulk.genes[x])
  #print(dim(current.genes))
  
  # get the numbers
  df <- data.frame(receptor = sum(pairs$receptor %in% current.genes),
                   ligand = sum(pairs$ligand %in% current.genes),
                   both = dim(pairs[pairs$receptor %in% current.genes & pairs$ligand %in% current.genes,])[1],
                   unique_receptor = length(unique(pairs[pairs$receptor %in% current.genes,1])),
                   unique_ligand = length(unique(pairs[pairs$ligand %in% current.genes,2])),
                   total_genes = length(current.genes)
  )
})

names(expression.per.celltype) <- celltypes
#View(expression.per.celltype)


# save results
numbers.report <- do.call(rbind.data.frame, expression.per.celltype)
#write.xlsx(numbers.report, file = paste(Sys.Date(), "Report_Q3.xlsx", sep = " "), sheetName = "Per_celltype", append = T)

#-----------------------------------------------------------------------------------------
### which receptors and ligands are present

interaction.matrix2 <- lapply(celltypes, function(x){
  # get current celltype
  current.celltype <- unlist(pseudobulk.genes[x])
  
  # check if gene is present from pair
  current.interactions <- sapply(rownames(pairs), function(i){
    # select the current receptor and current ligand and check presence
    receptor.expression <- pairs[i,"receptor"] %in% current.celltype
    ligand.expression <- pairs[i, "ligand"] %in% current.celltype
    
    # bind
    table <- cbind(receptor.expression, ligand.expression)
    names(table) <- c("receptor", "ligand")
    
    return(table)
    
  })
  
  return(current.interactions)
  
})

names(interaction.matrix2) <- celltypes

#-----------------------------------------------------------------------------------------
### interaction score heatmap
# get average expression per cluster of all possible ligands and receptors
receptors.and.ligands <- unique(c(as.character(pairs[,1]), as.character(pairs[,2])))

average.expression <-sapply(celltypes, function(type){
  cells <- WhichCells(seuset, ident = type)
  average <- gene.expression.matrix[receptors.and.ligands,cells]
  average <- apply(average, 1, mean)
})

# get average expression ligands

interaction.score <- lapply(celltypes, function(x){
  # get ligand info
  ligand.info <- data.frame(ligands.logical = interaction.matrix2[[x]]["ligand",],
                            ligands = as.character(pairs[,"ligand"]),
                            row.names = rownames(pairs), stringsAsFactors = F)
  #print(str(ligand.info))
  # get the ligand positions
  match.ligands <- match(ligand.info$ligands, rownames(average.expression))
  
  # add them to ligand info
  ligand.info$ligand.expression <- average.expression[match.ligands,x]
  #print(head(ligand.info))
  
  interaction <- lapply(celltypes, function(y){
    # get receptor info
    receptor.info <- data.frame(receptors.logical = interaction.matrix2[[y]]["receptor",],
                                receptors = as.character(pairs[,"receptor"]),
                                row.names = rownames(pairs), stringsAsFactors = F)
    #print(str(receptor.info))
    # get the receptor positions
    match.receptor <- match(receptor.info$receptors, rownames(average.expression))
    
    # add them to receptor info
    receptor.info$receptor.expression <- average.expression[match.receptor,y]
    #print(head(ligand.info))
    
    ## join both ligand and receptor
    res <- cbind(ligand.info, receptor.info)
    # calculate interaction score
    res$interaction.score <- res$ligand.expression * res$receptor.expression
    
    # # plot only T and T as defined earlier
    # plot <- res[res$ligands.logical ==T & res$receptors.logical ==T,]
    # #plot only interactions >= 1
    # #plot <- plot[ plot$interaction.score >= 1,]
    # # plot top 30 interactions
    # plot <- plot[order(plot$interaction.score, decreasing = T),]
    # 
    # 
    # #plot
    # pdf(paste(result.folder, "/", Sys.Date()," ", x," vs ", y, ".pdf", sep = ""))
    # pheatmap(plot[1:30,"interaction.score", drop = F], # drop = F to show rownames
    #          cluster_rows = F, # no clustering
    #          cluster_cols = F, 
    #          legend = T, 
    #          show_rownames = T, 
    #          show_colnames = F,
    #          main = paste(x,"vs", y, sep = " "))
    # dev.off()
    
    return(res)
  })
  
  names(interaction) <- celltypes
  return(interaction)
})

names(interaction.score) <- celltypes
View(interaction.score)

## plot interaction score distribution
pdf(paste(result.folder, "/", Sys.Date(), " interaction_score_histogram.pdf", sep = ""))
lapply(celltypes, function(x){
  current.list <- interaction.score[[x]]
  
  plot.it <- lapply(celltypes, function(y){
    # get only selected genes
    scores <- current.list[[y]]
    scores <- scores[scores$ligands.logical == T & scores$receptors.logical == T,]
    #head(scores)
    hist(scores$interaction.score, main = paste(x,y, sep = "_")) 
  })
})
dev.off()

#-----------------------------------------------------------------------------------------
### sum the interaction score

cumulative.interaction.score <- lapply(celltypes, function(x){
  # get current list
  current.list <- interaction.score[[x]]
  # get ligand
  ligand <- as.character(x)
  
  get.sum <- lapply(celltypes, function(y){
    # get current interaction table
    current.table <- current.list[[y]]
    
    # get ligand
    receptor <- as.character(y)
    
    # sum interaction score where ligands.logical and receptors.logical == T
    current.table <- current.table[current.table$ligands.logical == T & current.table$receptors.logical ==T,]
    cumulative.score <- sum(current.table$interaction.score)
    
    # write files
    #write.table(current.table[order(current.table$interaction.score, decreasing = T),], file = paste(ligand, receptor, "interaction.txt", sep = "_"))
    
    res <- cbind(ligand, receptor, cumulative.score)
    return(res)
  })
  
  res2 <- do.call(rbind.data.frame, get.sum)
  #print(str(res2))
  return(res2)
})



#View(cumulative.interaction.score)
cumulative.interaction.score <- do.call(rbind.data.frame, cumulative.interaction.score)
head(cumulative.interaction.score)
cumulative.interaction.score$cumulative.score <- as.numeric(as.character(cumulative.interaction.score$cumulative.score))

#-----------------------------------------------------------------------------------------
## plot cumulative scores celltype vs. celltype
# https://r.789695.n4.nabble.com/convert-data-frame-of-values-into-correlation-matrix-td1457592.html

plot.it <- cumulative.interaction.score
plot.it <- reshape2:: melt(plot.it)
plot.it <- reshape2:: dcast(plot.it, ligand ~ receptor)
rownames(plot.it) <- plot.it$ligand
plot.it <- plot.it[,2:24]
head(plot.it)
dim(plot.it)

pdf(paste(result.folder, "/", Sys.Date(), " cumulative_interaction_celltype_vs_celltype.pdf", sep = ""))
pheatmap::pheatmap(plot.it,
                   main = "Cumulative interaction scores")
dev.off()


### custom order
custom.order <- c("CD4.0", "CD4.1", "CD4.2", "CD4.3", "CD4.4", "CD8.0", "CD8.1", "CD8.2", "CD8.3", "E.0", "E.1", 
                  "E.2","E.3", "S.0", "S.1","My.0", "My.1", "My.2", "My.3", "My.4","KIT+ Mast Cells", "CD79A+ B-Cells", 
                  "Mixed Cells" )
plot.it2 <- plot.it[custom.order, custom.order]

numbers.report2 <- numbers.report[custom.order,]

pdf(paste(result.folder, "/", Sys.Date(), " pheatmap_cum_interaction_score.pdf", sep = ""))
# pheatmap rows = y axis = ligand
breaksList <- seq(0, max(plot.it2), by = 1)
pheatmap(plot.it2,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         main = "Cumulative interaction scores",
         breaks = breaksList,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white")
dev.off()


# # unique ligand receptor expression
# pdf(paste(result.folder, "/", Sys.Date(), " pheatmap_unique_ligand_receptor.pdf", sep = ""))
# pheatmap(numbers.report2[,c("unique_ligand", "unique_receptor")],
#          color =  rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Purples")))(100)),
#          cluster_rows = F,
#          cluster_cols = F,
#          border_color = "white")
# dev.off()

# ligand receptor expression
pdf(paste(result.folder, "/", Sys.Date(), " pheatmap_ligand_receptor.pdf", sep = ""))
breaksList <- seq(0, max(numbers.report2[,c("ligand", "receptor")]), by = 1)
pheatmap(numbers.report2[,c("ligand", "receptor")],
         color =  rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Purples")))(length(breaksList))),
         breaks = breaksList,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white")
dev.off()

pdf(paste(result.folder, "/", Sys.Date(), " pheatmap_ngene.pdf", sep = ""))
# total gene expression 
breaksList <- seq(0, max(numbers.report2[,"total_genes"]), by = 1)
pheatmap(numbers.report2[,"total_genes"],
         color =  rev(colorRampPalette(rev(brewer.pal(n = 7, name =  "Greens")))(length(breaksList))),
         breaks = breaksList,
         cluster_cols = F,
         cluster_rows = F,
         border_color = "white")
dev.off()


#-----------------------------------------------------------------------------------------
### Pathway analysis on top interactions
## top interactions only
# sort cumulative interaction score by size
top.cumulative.interaction.score <- cumulative.interaction.score[order(cumulative.interaction.score$cumulative.score, decreasing = T),]
top.cumulative.interaction.score$cellpairs <- paste(top.cumulative.interaction.score$ligand, 
                                                    top.cumulative.interaction.score$receptor,
                                                    sep = "_")

# worried that the 1:1 conversion from and to Entrez ID doesn't bring back the same results
# add column to pairs
Entrez.pairs <- data.frame(pairs,
                           Entrez.ligand = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(pairs$ligand), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID,
                           Entrez.receptor = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(pairs$receptor), columns = "ENTREZID", keytype = "SYMBOL", multiVals = "first")$ENTREZID )

# define the top x interactions
top.interactions <- 15

# get the genes and translate to entrezID
interacting.genes <- lapply(1:top.interactions, function(i){
  # get ligand expressing celltype
  from <- top.cumulative.interaction.score[i,"ligand"]
  # get receptor expressing celltype
  to <- top.cumulative.interaction.score[i,"receptor"]
  
  ## get the corresponding ligands
  ligands <- interaction.matrix2[[from]]["ligand",]
  ligands <- Entrez.pairs[ligands,"Entrez.ligand"]
  #ligands <- as.character(Entrez.pairs[ligands,"ligand"])
  
  ## get the corresponding receptors
  receptors <- interaction.matrix2[[to]]["receptor",]
  receptors <- Entrez.pairs[receptors, "Entrez.receptor"]
  #receptors <- as.character(Entrez.pairs[receptors, "receptor"])
  
  ## combine and unique per top interaction
  results <- unique(c(ligands,receptors))
  
  
  return(results)
})

names(interacting.genes) <- top.cumulative.interaction.score$cellpairs[1:top.interactions]
View(interacting.genes)

## pathway enrichment analysis
PEA <- lapply(1:top.interactions, function(x){
  # select genes
  current.genes <- interacting.genes[[x]]
  
  # pathway analysis 
  pathways <- enrichPathway(gene = current.genes, pvalueCutoff = 0.05, readable = T)
  
  #plot everyting
  name <- top.cumulative.interaction.score$cellpairs[x]
  
  #p1<- barplot(pathways, showCategory = 10)
  #p2<- dotplot(pathways, showCategory = 15)
  # ggsave(paste(result.folder, "/", Sys.Date(), "_", name, "_barpathway_analysis_.pdf", sep = ""), plot = p1)
  # ggsave(paste(result.folder, "/", Sys.Date(), "_", name,"_dotpathway_analysis_.pdf", sep = ""), plot = p2)
  # 
  
  p3 <- emapplot(pathways, showCategory = 25)
  
  ggsave(paste(result.folder, "/", Sys.Date(), "_", name, "_emaplot_.pdf", sep = ""), plot = p3)
  
  dev.off()
  
  return(as.data.frame(pathways))
})

View(PEA)

## picked result
cherry <- enrichPathway(gene = interacting.genes[["S.0_My2"]], pvalueCutoff = 0.05, readable = T)

mapplot(cherry, showCategory = 25)


#-----------------------------------------------------------------------------------------
### overlap with GWAS
# load file
enriched.GWAS <- read.table("enriched_GWAS.txt", header = T, stringsAsFactors = F)

# overlap with enriched?
sum(enriched.GWAS$Gene %in% receptors.and.ligands)

enriched.GWAS <- enriched.GWAS[enriched.GWAS$Gene %in% receptors.and.ligands,]
View(enriched.GWAS)
enriched.genes <- enriched.GWAS$Gene