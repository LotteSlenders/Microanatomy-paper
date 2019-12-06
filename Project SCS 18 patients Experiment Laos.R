### 2019-5-16

# Project SCS 18 patients Experiment Laos

# single cell vs Bulk

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Laos")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(plyr)          # v0.7.8
library(SingleR)        # v0.2.2
library(org.Hs.eg.db)

# object
load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined
rm(all.seur.combined)

# set seed
set.seed(69)

#-----------------------------------------------------------------------------------------
### create count matrix for single cells

## prep meta data and files
# load meta data
meta.data <- read.table('Project SCS 18 patients meta data.txt', header = T, skipNul = T, fill = T, stringsAsFactors = F)
# rename first column
names(meta.data)[1] <- "Patient"
# add ID column
meta.data$ID <- paste(meta.data$Patient, ".P",meta.data$Plate, sep = "")
# transform patient colum from integer to character
meta.data$Patient <- as.character(meta.data$Patient)

sum(duplicated(meta.data$ID))
#View(meta.data)

# check files
all.files <- list.files(pattern="*.tsv")
all.files
all.files %in% meta.data$File
length(all.files) == dim(meta.data)[1]

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

library(org.Hs.eg.db)   # v3.7.0

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


# merge suerat objects from list
seuset <- seurat.object.list[[1]]
for (i in 2:length(x = seurat.object.list)) {
    seuset <- MergeSeurat(object1 = seuset, object2 = seurat.object.list[[i]], do.normalize = F)
}

dim(seuset@data)

single.cell <- as.matrix(seuset@data)
single.cell <- data.frame(Pseudobulk = rowSums(single.cell), row.names = rownames(single.cell))
# write.csv(single.cell, "Single_cell_all_genes.txt")

rm(seurat.object.list)

#-----------------------------------------------------------------------------------------
### create count matrix for bulk
whole.tissue <- read.delim("DH_combined_counts.names.txt", header = TRUE)
dim(whole.tissue)
head(whole.tissue)
# remove unwanted genes
unwanted.genes <- c("^ERCC", "^MT-", "^RPL", "^RPS", "^KCNQ1OT1$", "^UGDH-AS1$", "^MALAT1$", "^EEF1A$")
whole.tissue <- whole.tissue[grep(paste(unwanted.genes, collapse = "|"), whole.tissue$X ,invert=TRUE),]
dim(whole.tissue)

# update gene symbols
update.symbols <- function(my.counts.table = my.counts.table) {
  #Keytype has to be alias to retrieve the 'old' symbols, otherwise will get NAs instead of updated names.
  #Retrieve the mapping
  current.keys <- as.character(my.counts.table$X)
  gene.info <- AnnotationDbi::select(org.Hs.eg.db, keys = current.keys, columns = "SYMBOL", keytype = "ALIAS", multiVals = "first")
  
  #Keep only one mapping per alias
  d <- duplicated(gene.info$ALIAS)
  gene.info <- gene.info[!d,]
  
  #Add mappings to the table
  my.counts.table$SYMBOL <- my.counts.table$X
  my.counts.table$SYMBOL <- mapvalues(my.counts.table$SYMBOL, from=gene.info$ALIAS, to=gene.info$SYMBOL)
  
  
  #Remove non-mappings (old LOCs and stuff that are not availble anymore)
  na <- is.na(my.counts.table$SYMBOL)
  my.counts.table <- my.counts.table[!na,]
  
  # pick a random duplicated entry 
  my.counts.table <- my.counts.table %>% distinct(SYMBOL, .keep_all = TRUE) 
  
}

# update gene names
whole.tissue <- update.symbols(whole.tissue)
# set new rownames
rownames(whole.tissue) <- whole.tissue$SYMBOL
head(whole.tissue)

# select colums
whole.tissue <- whole.tissue[,2:3]
names(whole.tissue) <- c("Bulk1", "Bulk2")
head(whole.tissue)

#-----------------------------------------------------------------------------------------
### join two datasets

# add common variables for joining
single.cell$SYMBOL <- rownames(single.cell)
whole.tissue$SYMBOL <- rownames(whole.tissue)
# join
single.vs.whole <- full_join(single.cell, whole.tissue, by = "SYMBOL")
# clean up
rownames(single.vs.whole) <- single.vs.whole$SYMBOL
head(single.vs.whole)
single.vs.whole <- single.vs.whole[,c("Pseudobulk", "Bulk1", "Bulk2")]
single.vs.whole <- round(single.vs.whole)
single.vs.whole[is.na(single.vs.whole)] <- 0
head(single.vs.whole)

write.table(single.vs.whole, "single_vs_whole.txt")
#-----------------------------------------------------------------------------------------
### get top 10 markers
## get marker genes/top differentially expressed
# note that marker genes are only positively DEG
#load and save
all.marker.files <- list.files(pattern="*genes.txt")
all.marker.files
marker.genes.list <- list()

for (i in seq_along(all.marker.files)) {
  current.markers <- read.delim(all.marker.files[i], stringsAsFactors = F)
  current.markers <- current.markers[1:10,"gene"]
  marker.genes.list[[i]] <- current.markers
}

names(marker.genes.list) <- paste("Cluster", 0:13, sep = " ")
View(marker.genes.list)


#------------------------------------------------------------------------------
### plot
## pseudobulk vs bulk1
# get top 10's
all.top.10 <- reshape2::melt(marker.genes.list)
all.top.10 <- unlist(all.top.10$value)

# define plotting data
plot.pseudo <- log2(single.vs.whole[,c(1,2)] + 0.1)

# get colours
colour.palette <- scales::hue_pal()(14)

#plot
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2]), 
             colour = "red", size = 4) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("Top 10 DEG of each celltype", sep = " ")) 
  #geom_text(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2], label = all.top.10))
ggsave(paste(result.folder, "/", Sys.Date(), "Pseudobulk vs. bulk1 all marker genes in red.pdf", sep = ""))

## all marker genes

lapply(marker.genes.list, function(x){
  ggplot(data = plot.pseudo, aes(x = Pseudobulk, y = Bulk1)) +
  geom_point(colour = "black", size = 1) +
  geom_point(data = plot.pseudo[x,], inherit.aes = T, colour = "red", size = 4) +
    ggtitle()
  
})

lapply(marker.genes.list, function(x) {
  print(plot.pseudo[x,])
})

single.cell["PGRM5P2",]

#"PGM5P2", "EEF1A1", "MAB21L3"

# get the top 10's
top.10s <- as.data.frame(marker.genes.list, stringsAsFactors = F)
top.10s[,2]
names(top.10s) <- paste("Cluster", 0:13, sep = " ")
head(top.10s)


ggplot(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2])) + 
  geom_point( colour = "black", size = 1) +
  # CD3+CD8+ T-Cells I
  geom_point(data = plot.pseudo[top.10s[,1],], aes(plot.pseudo[top.10s[,1],1], plot.pseudo[top.10s[,1],2], label = top.10s[,1]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[1]) + geom_text()
  # CD3+CD4+ T-Cells I
  geom_point(data = plot.pseudo[top.10s[,2],], aes(plot.pseudo[top.10s[,2],1], plot.pseudo[top.10s[,2],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[2]) +
  # Mixed cells
  geom_point(data = plot.pseudo[top.10s[,6],], aes(plot.pseudo[top.10s[,6],1], plot.pseudo[top.10s[,6],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[6] ) +
  # CD3+CD8+ T-Cells II
  geom_point(data = plot.pseudo[top.10s[,4],], aes(plot.pseudo[top.10s[,4],1], plot.pseudo[top.10s[,4],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[4]) +
  # ACD3+CD4+ T-Cells II
  geom_point(data = plot.pseudo[top.10s[,3],], aes(plot.pseudo[top.10s[,3],1], plot.pseudo[top.10s[,3],2]),
             colour = "black", shape = 21, fill = colour.palette[3], size = 5) +
  # CD14+CD68+ Macrophages I
  geom_point(data = plot.pseudo[top.10s[,5],], aes(plot.pseudo[top.10s[,5],1], plot.pseudo[top.10s[,5],2]),
             colour = "black", shape = 21, fill = colour.palette[5], size = 5) +
  # CD14+CD68+ Macrophages II
  geom_point(data = plot.pseudo[top.10s[,7],], aes(plot.pseudo[top.10s[,7],1], plot.pseudo[top.10s[,7],2]),
             colour = "black", shape = 21, fill = colour.palette[7], size = 5) +
  # CD14+CD68+ Macrophages III
  geom_point(data = plot.pseudo[top.10s[,8],], aes(plot.pseudo[top.10s[,8],1], plot.pseudo[top.10s[,8],2]),
             colour = "black", shape = 21, fill = colour.palette[8], size = 5) +
  # MYH11+ Smooth muscle cells
  geom_point(data = plot.pseudo[top.10s[,9],], aes(plot.pseudo[top.10s[,9],1], plot.pseudo[top.10s[,9],2]),
             colour = "black", shape = 21, fill = colour.palette[9], size = 5) +
  # CD34+ endothelial Cells I
  geom_point(data = plot.pseudo[top.10s[,10],], aes(plot.pseudo[top.10s[,10],1], plot.pseudo[top.10s[,10],2]),
             colour = "black", shape = 21, fill = colour.palette[10], size = 5) +
  # CD34+ endothelial Cells II
  geom_point(data = plot.pseudo[top.10s[,11],], aes(plot.pseudo[top.10s[,11],1], plot.pseudo[top.10s[,11],2]),
             colour = "black", shape = 21, fill = colour.palette[11], size = 5) +
  # CD79A+ B cells
  geom_point(data = plot.pseudo[top.10s[,12],], aes(plot.pseudo[top.10s[,12],1], plot.pseudo[top.10s[,12],2]),
             colour = "black", shape = 21, fill = colour.palette[12], size = 5) +
  # CD14+CD68+ Macrophages IV
  geom_point(data = plot.pseudo[top.10s[,13],], aes(plot.pseudo[top.10s[,13],1], plot.pseudo[top.10s[,13],2]),
             colour = "black", shape = 21, fill = colour.palette[13], size = 5) +
  # KIT+ Mast Cells
  geom_point(data = plot.pseudo[top.10s[,14],], aes(plot.pseudo[top.10s[,14],1], plot.pseudo[top.10s[,14],2]),
             colour = "black", shape = 21, fill = colour.palette[14], size = 5) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2]))


ggplot(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2])) + 
  geom_point( colour = "black", size = 1) +
  geom_point(data = plot.pseudo[top.10s[3,2],], aes(plot.pseudo[top.10s[3,2],1], plot.pseudo[top.10s[3,2],2]),
             colour = "black", size = 5, shape = 21, fill = colour.palette[2])

single.vs.whole["ADGRG1",]  
plot.pseudo["ADGRG1",]
