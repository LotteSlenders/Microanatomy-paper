### 2019-1-23

# Project Single Cell Sequencing Experiment Heartstrings

# Figure and script for Single cell paper Figure 1

#-----------------------------------------------------------------------------------------
# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project Single Cell Plaque/Experiment Heartstrings")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.1

### load packages
library(Seurat)         # V2.3.4
library(dplyr)          # v0.7.6
library(circlize)       # v0.4.5
library(org.Hs.eg.db)
library(AnnotationDbi)
library(plyr)
library(tidyr)
library(DESeq2)

#-----------------------------------------------------------------------------------------
### Expressed in bulk
# Expression differences between whole tissue and singe cell

### Prep single cell data
#load
all.files <- list.files(pattern="*.tsv")
single.cell <- read.delim(all.files[1])
for (df in c(2:length(all.files))) {
  next.df <- read.delim(all.files[df])
  single.cell <- full_join(single.cell, next.df, by = "GENEID")
}

# replace NA with 0
single.cell[is.na(single.cell)] <- 0
# fix GENEID (get rid of chromosome info)
single.cell$GENEID <- sub("__.*", "" , single.cell$GENEID)
# Make it pretty
rownames(single.cell) <- single.cell$GENEID
single.cell <- single.cell[,-1]

## sum genecounts per patient
number.of.plates.per.patient <- 3
number.of.plates <- length(all.files)
number.of.patients <- number.of.plates / number.of.plates.per.patient
cells.per.plate <- 384
single.cell.sum <- data.frame(row.names = rownames(single.cell))


for (patient in c(1:number.of.patients)) {
  # create variable named after patient
  #new.name <- paste("sum.P", patient, sep = "")
  # subset data per patient
  plates.end <- (patient * number.of.plates.per.patient) * cells.per.plate
  plates.start <- plates.end - (number.of.plates.per.patient * cells.per.plate) + 1
  # sum counts
  sum.per.patient <- rowSums(single.cell[,plates.start:plates.end])
  #add to dataframe
  single.cell.sum[,patient] <- sum.per.patient
  
}

# calculate the average
single.cell.sum$average <- rowMeans(single.cell.sum)
# calculate a pseudobulk
single.cell.sum$sum <- rowSums(single.cell.sum[,1:3])

# add common variable for joining
single.cell.sum$X <- rownames(single.cell)
dim(single.cell.sum)

### Prep whole tissue data
whole.tissue <- read.delim("DH_combined_counts.names.txt", header = TRUE)
dim(whole.tissue)

### merge single cell and whole tissue
single.vs.whole <- full_join(whole.tissue[,c(2:4)], single.cell.sum, by = "X")
colnames(single.vs.whole) <- c("Bulk1", "Bulk2", "GENEID", "Single1", "Single2", "Single3", "Average single", "Pseudobulk")
dim(single.vs.whole)
head(single.vs.whole)

# remove unwanted genes and genes without gene ID
unwanted.genes <- c("ERCC", "MT-","RPL", "RPS","UGDH-AS1","KCNQ1OT1", "MALAT1", "EEF1A1")
single.vs.whole <- single.vs.whole[grep(paste(unwanted.genes, collapse = "|"),single.vs.whole$GENEID ,invert=TRUE),]
dim(single.vs.whole)
# remove rows with empty GENEID
single.vs.whole  <-single.vs.whole[single.vs.whole$GENEID != "",]
dim(single.vs.whole)
# replace NA with 0
single.vs.whole[is.na(single.vs.whole)] <- 0

# remove al duplicate entries
single.vs.whole <- single.vs.whole[!duplicated(single.vs.whole$GENEID),]
dim(single.vs.whole)
head(single.vs.whole)

# now there are no more duplicates and GENEID an be set as rownames
rownames(single.vs.whole) <- single.vs.whole$GENEID
# Make it pretty
single.vs.whole <- single.vs.whole[,c("Bulk1", "Bulk2", "Single1", "Single2", "Single3", "Average single", "Pseudobulk")]
dim(single.vs.whole)
head(single.vs.whole)

plot(log2(single.vs.whole[,c(1,7)]+1), pch = 1, cex = 0.5)
cor(single.vs.whole)
# > plot(log2(single.vs.whole[,c(2,1)]+1))
cor(log2(single.vs.whole+ 0.1))

pdf(paste(result.folder, "/", Sys.Date(), " Bulk vs single cells.pdf", sep = ""), width = 10)
plot(log2(single.vs.whole + 0.1), pch = 1, cex = 0.07)
dev.off()


### remove genes that are not expressed in one of the two datasets
# prepare all possible combinations
choices <- choose(7, 2)       # identify the number of possible pairs
combinations <- combn(7, 2)    # create those pairs
correlations <- c()

for (choice in 1:choices) {
  # select the cluster numbers to run
  dataset1 <- combinations[,choice][1] # the first index of the combination
  dataset2 <- combinations[,choice][2] # the second index of the combination
  
  # subset data
  df <- single.vs.whole[,c(dataset1, dataset2)]
  
  # save plot
  pdf(paste(result.folder, "/", Sys.Date(), " ", colnames(single.vs.whole[combinations[,choice][1]]), " vs.", 
            colnames(single.vs.whole[combinations[,choice][2]]), ".pdf", sep = ""))
  plot(log2(df + 0.1), pch = 1, cex = 0.07)
  dev.off()
  
  # remove all instances where either dataset1 or dataset 2 has no entries
  df <- df[df[,1] != 0 & df[,2] != 0,]
  
  # correlation
  current.correlation <- cor(df[,1], df[,2])
  correlations <- c(correlations, current.correlation)
}

#-----------------------------------------------------------------------------------------
### Prep whole tissue data
whole.tissue <- read.delim("DH_combined_counts.names.txt", header = TRUE)
dim(whole.tissue)
# remove unwanted genes
unwanted.genes <- c("^ERCC", "^MT-", "^RPL", "^RPS", "^KCNQ1OT1$", "^UGDH-AS1$", "^MALAT1$", "^EEF1A$")
whole.tissue <- whole.tissue[grep(paste(unwanted.genes, collapse = "|"), whole.tissue$X ,invert=TRUE),]
dim(whole.tissue)

# split in Bulk1 and Bulk2
Bulk1 <- whole.tissue[,c(1,2,4)]
Bulk2 <- whole.tissue[,c(1,3,4)]
# remove gene entries with 0 counts
Bulk1 <- Bulk1[Bulk1[,2] !=0,]
dim(Bulk1)
Bulk2 <- Bulk2[Bulk2[,2] != 0,]
dim(Bulk2)

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
    
  # pick a random duplicated entry since I only really need the gene names
  my.counts.table <- my.counts.table %>% distinct(SYMBOL, .keep_all = TRUE) 
   
}

Bulk1 <- update.symbols(Bulk1)
Bulk2 <- update.symbols(Bulk2)
rownames(Bulk1) <- Bulk1$SYMBOL

# single cell data
seuset <- readRDS("all.seur.combined.RDS")

## meta information
# Fetch cell identity
ID <- FetchData(object = seuset,
                vars.all = c("ident", "orig.ident"))
dim(ID)
# Fetch meta data
meta.data <- seuset@meta.data
dim(meta.data)

# combine 
meta.data <- cbind(meta.data, ID)
# check
dim(meta.data)
sum(rownames(meta.data) == meta.data$orig.ident)

# count data
# count matrix, normalized and log-transformed single cell expression.
# for data without preprocessing use object@raw.data
counts <- as.matrix(seuset@raw.data)
dim(counts)
# check
sum(rownames(meta.data) == colnames(counts))

# celltypes
celltypes <- levels(unique(meta.data$ident))
celltypes

### subset per celltype
counts.per.celltype <- list()

for (i in seq_along(celltypes)) {
  # get celltype
  celltype <- celltypes[i]
  # select only wells corresponding to celltype
  wells <- meta.data[meta.data$ident == celltype, "orig.ident"]
  # select only counts corresponding to celltype
  current.counts <- counts[,wells]
  # save counts per celltype
  counts.per.celltype[[i]] <- current.counts

}

### Circosplot
my.circos <- function(names, bulk1.in.type, type.in.bulk1, colors, bg.colors){
  # Ensure that each sector has a domain of 0 to 1 (x) and a non-zero range (y)
  # every sector has its own name, given by 'names'
  data = data.frame(
    factor = names,
    x = c(rep(0, length(names)), rep(1, length(names))),
    y = c(0,1)
  ) 
  
  # Initialize the plot. (mar = margin)
  par(mar = c(1, 1, 1, 1) ) 
  circos.initialize(factors = data$factor, x = data$x )
  
  # Build the regions of track #1
  #circos.trackPlotRegion(factors = data$factor, y=data$y , bg.col = bg.colors , bg.border = NA )
                         #,circos.axis(labels = T))
  
  circos.trackPlotRegion(factors = data$factor, y=data$y , bg.col = bg.colors , bg.border = NA, panel.fun = function(x,y) {
    circos.text(x = 0.5, y = 1.1,
      labels = names,
      facing = "outside",
      niceFacing = T
      
    )
  } )
  
  # For each cell type
  for(i in seq(2, length(names))){
    # Get identity
    name <- names[i] 
    # We want to draw the area from the middle of the sectors.
    # The middle is at 0.5, hence we + and - the half of the percentage from this middle to center
    # the link nicely in the sector.
    percentage.from <- bulk1.in.type[i] / 2
    percentage.to <- type.in.bulk1[i] / 2
    color <- colors[i]
    circos.link(name, c(0.5 - percentage.from, 0.5 + percentage.from), "bulk1", 
                c(0.5 - percentage.to, 0.5 + percentage.to),
                col = color)
    
  } 
}


names <- c('bulk1')
bulk1.in.type <- c(1)
type.in.bulk1 <- c(1)
bulk.rownames = rownames(Bulk1)
for(i in seq_along(celltypes)){
  cell.type <- counts.per.celltype[[i]]
  cleansed <- cell.type[rowSums(cell.type)>0, ]
  type.rownames <- rownames(cleansed) 
  this.type.in.bulk <- sum(type.rownames %in% bulk.rownames) 
  this.bulk.in.type <- sum(bulk.rownames %in% type.rownames) 
  type.in.bulk1 <- c(type.in.bulk1, (this.type.in.bulk) / length(bulk.rownames))
  bulk1.in.type <- c(bulk1.in.type, (this.bulk.in.type) / length(type.rownames)) 
  names <- c(names, celltypes[i])
}

#library("RColorBrewer")
#colors <- brewer.pal(n = length(names), name = 'RdBu')
colors <- rainbow(length(names), alpha = 0.3)
bg.colors <- rainbow(length(names), alpha = 0.7)
#colors <- rep(rgb(0, 1, 0, alpha=0.3), length(names))
my.circos(names, bulk1.in.type, type.in.bulk1, colors, bg.colors)



name=c(3,10,10,3,6,7,8,3,6,1,2,2,6,10,2,3,3,10,4,5,9,10)
feature=paste("feature ", c(1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5) , sep="")
dat <- data.frame(name,feature)
dat <- with(dat, table(name, feature))

# Charge the circlize library
library(circlize)

# Make the circular plot
chordDiagram(as.data.frame(dat), transparency = 0.5)


names <- c('bulk1', 't1', 't2', 't3')
bulk1.in.type <- c(1, 0.6, 0.2, 0.75)
type.in.bulk1 <- c(1, 0.8, 0.5, 0.4)
colors <- rep(rgb(0, 1, 0, alpha=0.3), length(names))

#-----------------------------------------------------------------------------------------
### try again
### Prep whole tissue data
whole.tissue <- read.delim("DH_combined_counts.names.txt", header = TRUE)
dim(whole.tissue)
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

whole.tissue <- update.symbols(whole.tissue)
rownames(whole.tissue) <- whole.tissue$X

### prep single cell data
seuset <- readRDS("all.seur.combined.RDS")

# Fetch cell identity
ID <- FetchData(object = seuset,
                vars.all = c("ident", "orig.ident"))
dim(ID)
# Fetch meta data
meta.data <- seuset@meta.data
dim(meta.data)

# combine 
meta.data <- cbind(meta.data, ID)
# check
dim(meta.data)
sum(rownames(meta.data) == meta.data$orig.ident)

# count data
# count matrix, normalized and log-transformed single cell expression.
# for data without preprocessing use object@raw.data
counts <- as.matrix(seuset@raw.data)
dim(counts)
# check
sum(rownames(meta.data) == colnames(counts))

# celltypes
celltypes <- levels(unique(meta.data$ident))
celltypes

## subset per celltype
counts.per.celltype <- list()
sum.per.celltype <- data.frame(matrix(ncol = length(celltypes), nrow = length(rownames(seuset@raw.data))))
rownames(sum.per.celltype) <- rownames(seuset@raw.data)

for (i in seq_along(celltypes)) {
  # get celltype
  celltype <- celltypes[i]
  # select only wells corresponding to celltype
  wells <- meta.data[meta.data$ident == celltype, "orig.ident"]
  # select only counts corresponding to celltype
  current.counts <- counts[,wells]
  # save counts per celltype
  counts.per.celltype[[i]] <- current.counts
  
  # save pseudobulk per celltype
  sums <- rowSums(current.counts)
  sum.per.celltype[,i] <- sums
  
}


## combine single and bulk
# add common value for joining
sum.per.celltype$SYMBOL <- rownames(sum.per.celltype)
# join by SYMBOL
single.vs.whole <- full_join(sum.per.celltype, whole.tissue[,c(2,3,5)], by = "SYMBOL")
head(single.vs.whole)
# set rownames again
rownames(single.vs.whole) <- single.vs.whole$SYMBOL
# set colnames
names(single.vs.whole) <- c(celltypes, "SYMBOL", "Bulk1", "Bulk2")
head(single.vs.whole)
# add pseudobulk
single.vs.whole$pseudobulk <- rowSums(single.vs.whole[,1:length(celltypes)])
# remove SYMBOL
single.vs.whole <- single.vs.whole[,c(1:11,13,14,15)]
head(single.vs.whole)
# replace NA with 0
single.vs.whole[is.na(single.vs.whole)] <- 0

#write.table(single.vs.whole, "Single cell vs Bulk counts.txt", sep = ",")

### single vs bulk
## plot bulk vs single cell and highlight top x marker genes per celltype

#plot(log2(single.vs.whole + 0.1), pch = 1, cex = 0.07)
cor(log2(single.vs.whole + 0.1))
#write.table(cor(log2(single.vs.whole + 0.1)), file = "correlations.txt", sep = ",", row.names = T, col.names = T)

## get marker genes/top differentially expressed
# note that marker genes are only positively DEG
#load and save
all.marker.files <- list.files(pattern="*genes.txt")
marker.genes.list <- list()

for (i in seq_along(all.marker.files)) {
  current.markers <- read.delim(all.marker.files[i])
  marker.genes.list[[i]] <- current.markers
}

# get same colours as used in seurat (https://github.com/satijalab/seurat/issues/257)
colour.palette <- scales::hue_pal()(length(celltypes))

### plot
# loading of marker files is not sequential
file.order <- readr::parse_number(all.marker.files) + 1
all.top.10 <- c()
top.10s <- data.frame(matrix(nrow = 10))

for (i in 1:length(celltypes)) {
  # for each cellype, get the correct marker genes
  # get the correct list by seeing which file corresponds to the celltype
  list.number <- which(file.order == i)
  current.df <- data.frame(marker.genes.list[[list.number]])
  # get topgenes
  top.genes <- levels(current.df$gene)
  # get top 10 (it only works with a two step) and store
  top.10 <- top.genes[1:10]
  all.top.10 <- c(all.top.10, top.10)
  
  # top 10 list
  top.10s[,i] <- top.10
  
  ## plot and save
  # plot i vs bulk1
  # plot.bulk <- log2(single.vs.whole[,c(12,i)] + 0.1)
  # ggplot() +
  #   geom_point(data = plot.bulk, aes(plot.bulk[,1], plot.bulk[,2]), colour = "black", size = 1) +
  #   geom_point(data = plot.bulk[top.10,], aes(plot.bulk[top.10,1], plot.bulk[top.10,2]), colour = colour.palette[i], size = 4) +
  #   xlab(names(plot.bulk[1])) + ylab(names(plot.bulk[2])) +
  #   ggtitle(paste("Top 10 DEG", celltypes[i], sep = " "))
  # #   geom_text(data = plot.bulk[top.10,], aes(plot.bulk[top.10,1], plot.bulk[top.10,2]), label = top.10, col = "red", hjust = 0, vjust = 0)
  # ggsave((paste(result.folder, "/", Sys.Date(), " ", celltypes[i], " vs. bulk1.pdf", sep = "")))
  # 
  # # plot i vs bulk 2
  # plot.bulk <- log2(single.vs.whole[,c(13,i)] + 0.1)
  # ggplot() +
  #   geom_point(data = plot.bulk, aes(plot.bulk[,1], plot.bulk[,2]), colour = "black", size = 1) +
  #   geom_point(data = plot.bulk[top.10,], aes(plot.bulk[top.10,1], plot.bulk[top.10,2]), colour = colour.palette[i], size = 4) +
  #   xlab(names(plot.bulk[1])) + ylab(names(plot.bulk[2])) +
  #   ggtitle(paste("Top 10 DEG", celltypes[i], sep = " "))
  # #   geom_text(data = plot.bulk[top.10,], aes(plot.bulk[top.10,1], plot.bulk[top.10,2]), label = top.10, col = "red", hjust = 0, vjust = 0)
  # ggsave(paste(result.folder, "/", Sys.Date(), " ", celltypes[i], " vs. bulk2.pdf", sep = ""))
  
  
  # top 10 in speudobulk vs bulk1
  # plot.bulk <- log2(single.vs.whole[,c(12,14)] + 0.1)
  # ggplot(data = plot.bulk, aes(plot.bulk[,1], plot.bulk[,2])) +
  #   geom_point(colour = "black", size = 1) +
  #   geom_point(data = plot.bulk[top.10,], aes(plot.bulk[top.10,1], plot.bulk[top.10,2]), colour = colour.palette[i], size = 4) +
  #   xlab("Bulk") + ylab(names(plot.bulk[2])) +
  #   ggtitle(paste("Top 10 DEG", celltypes[i], sep = " ")) 
  # #+ geom_smooth(method = 'lm', formula=y~x-1, se = F, col = "red") +
  # ggsave((paste(result.folder, "/", Sys.Date()," ",i, " PSeudobulk vs. bulk1 ", celltypes[i]," .pdf", sep = "")))
  
}

## plot pseudobulk vs bulk with all top 10
# pseudobulk vs bulk1
plot.pseudo <- log2(single.vs.whole[,c(12,14)] + 0.1)
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2]), 
             colour = "red", size = 4) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("Top 10 DEG of each celltype", sep = " "))  
  #+ geom_text(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2], label = all.top.10))
ggsave(paste(result.folder, "/", Sys.Date(), "Pseudobulk vs. bulk1.pdf", sep = ""))

# pseudobulk vs bulk2
plot.pseudo <- log2(single.vs.whole[,c(13,14)] + 0.1)
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2]), 
             colour = "red", size = 4) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("Top 10 DEG of each celltype", sep = " "))  
ggsave(paste(result.folder, "/", Sys.Date(), "Pseudobulk vs. bulk2.pdf", sep = ""))


### get genes present in SC but not Bulk and vise versa
not.expressed <- rownames(single.vs.whole[single.vs.whole$pseudobulk == 0,])
#write.table(not.expressed, "not expressed in SC data.txt")

is.expressed <- rownames(single.vs.whole[single.vs.whole$pseudobulk != 0 & single.vs.whole$Bulk1 == 0 | single.vs.whole$Bulk2 == 0,])
#write.table(is.expressed, "expressed in SC but not in bulk.txt")

### DEG
coldata <- data.frame(name = names(single.vs.whole[,c(12:14)]),
                      type = c("Bulk", "Bulk", "Pseudobulk"))
dds <- DESeqDataSetFromMatrix(countData =  single.vs.whole[,c(12:14)],
                               colData = coldata,
                               design = ~ type)
dds <- DESeq(dds)                               
pre.results <- results(dds)
head(pre.results)
plotMA(pre.results)

## select DEG
dim(pre.results)
results <- pre.results[!is.na(pre.results$padj),]
dim(results)
results <- results[results$padj <= 0.05,]
dim(results)
# sort by p value
results <- results[order(results$padj),]
# get top x positive and negative
top.pos.DEG <- rownames(results[results$log2FoldChange >= 0,])[1:100]
top.neg.DEG <- rownames(results[results$log2FoldChange <= 0,])[1:100]

# save
#write.table(results[results$log2FoldChange <= 0,], "significant negative DEG.txt")
#write.table(results[results$log2FoldChange >= 0,], "significant positive DEG.txt")

# plot
## plot pseudobulk vs bulk top DEG
# pseudobulk vs bulk1
plot.pseudo <- log2(single.vs.whole[,c(12,14)] + 0.1)
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[top.pos.DEG,], aes(plot.pseudo[top.pos.DEG,1], plot.pseudo[top.pos.DEG,2]), 
             colour = "red", size = 2) +
  geom_point(data = plot.pseudo[top.neg.DEG,], aes(plot.pseudo[top.neg.DEG,1], plot.pseudo[top.neg.DEG,2]), 
             colour = "blue", size = 2) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("200 DEG pseudobulk vs Bulk", sep = " "))  
ggsave(paste(result.folder, "/", Sys.Date(), " 200 DEG Pseudobulk vs. bulk1.pdf", sep = ""))

# pseudobulk vs bulk2
plot.pseudo <- log2(single.vs.whole[,c(13,14)] + 0.1)
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[top.pos.DEG,], aes(plot.pseudo[top.pos.DEG,1], plot.pseudo[top.pos.DEG,2]), 
             colour = "red", size = 2) +
  geom_point(data = plot.pseudo[top.neg.DEG,], aes(plot.pseudo[top.neg.DEG,1], plot.pseudo[top.neg.DEG,2]), 
             colour = "blue", size = 2) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("200 DEG pseudobulk vs Bulk", sep = " "))  
ggsave(paste(result.folder, "/", Sys.Date(), " 200 DEG Pseudobulk vs. bulk2.pdf", sep = ""))

# plot with top 10 all behind it
plot.pseudo <- log2(single.vs.whole[,c(12,14)] + 0.1)
ggplot() + 
  geom_point(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2] ), colour = "black", size = 1) +
  geom_point(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2]), 
            colour = "green", size = 3) +
  geom_point(data = plot.pseudo[top.pos.DEG,], aes(plot.pseudo[top.pos.DEG,1], plot.pseudo[top.pos.DEG,2]), 
             colour = "red", size = 2) +
  geom_point(data = plot.pseudo[top.neg.DEG,], aes(plot.pseudo[top.neg.DEG,1], plot.pseudo[top.neg.DEG,2]), 
             colour = "blue", size = 2) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("200 DEG pseudobulk vs Bulk and all marker genes", sep = " "))  
ggsave(paste(result.folder, "/", Sys.Date(), " 200 DEG + marker genes Pseudobulk vs. bulk1.pdf", sep = ""))

# plot whatever
plot.pseudo <- log2(single.vs.whole[,c(12,14)] + 0.1)
input <- rownames(results[results$log2FoldChange >= 0,])
ggplot(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2])) + 
  geom_point(colour = "black", size = 1) +
  geom_point(data = plot.pseudo[all.top.10,], aes(plot.pseudo[all.top.10,1], plot.pseudo[all.top.10,2]), 
             colour = "green", size = 3) +
  geom_point(data = plot.pseudo[input,], aes(plot.pseudo[input,1], plot.pseudo[input,2]), 
             colour = "red", size = 2) +
  xlab(names(plot.pseudo[1])) + ylab(names(plot.pseudo[2])) +
  ggtitle(paste("significant DEG positive foldchange (padj <0.05) + marker genes", sep = " "))  
ggsave(paste(result.folder, "/", Sys.Date(), " marker genes + significant DEG positive foldchange Pseudobulk vs. bulk1.pdf", sep = ""))


#-----------------------------------
### estimate population of T cell genes in Bulk based on marker genes
## compare t cells specific gene expression between singlecell and bulk
## https://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists

# overrepresentation function
get.overrepresentation <- function(marker.genes) {
  # prep data                  
  pseudobulk.marker.sum <- sum(single.vs.whole[marker.genes,14])
  Bulk1.marker.sum <- sum(single.vs.whole[marker.genes,12])
  Bulk2.marker.sum <- sum(single.vs.whole[marker.genes,13])
  pseudobulk.rest <- sum(single.vs.whole$pseudobulk) - pseudobulk.marker.sum
  Bulk1.rest <- sum(single.vs.whole$Bulk1) - Bulk1.marker.sum
  Bulk2.rest <- sum(single.vs.whole$Bulk2) - Bulk2.marker.sum
  
  # run statistical test
  print(fisher.test(matrix(c(pseudobulk.marker.sum,pseudobulk.rest,Bulk1.marker.sum,Bulk1.rest),nrow=2,ncol=2),alternative="greater"))
  print(fisher.test(matrix(c(pseudobulk.marker.sum,pseudobulk.rest,Bulk2.marker.sum,Bulk2.rest),nrow=2,ncol=2),alternative="greater"))
  
  # get folddifference
  prop1 <- pseudobulk.marker.sum / sum(single.vs.whole$pseudobulk)
  prop2 <- Bulk1.marker.sum / sum(single.vs.whole$Bulk1)
  prop3 <- Bulk2.marker.sum / sum(single.vs.whole$Bulk2)
  
  print(prop1/prop2)
  print(prop1/prop3)                    
                    
}

# get markers function
get.markers <- function(file.numbers) {
  current.markers <- c()
  for (i in file.numbers) {
  current.marker.genes <- data.frame(marker.genes.list[[i]])
  current.marker.genes <- current.marker.genes[current.marker.genes$p_val_adj <= 0.05,]
  current.marker.genes <- levels(current.marker.genes$gene)
  current.markers <- c(current.markers, current.marker.genes)
  }
  return(unique(current.markers))
} 

# t cells
tcell.marker.genes <- get.markers(c(1,2,5,7))
get.overrepresentation(tcell.marker.genes)

# macrophages
macrophage.marker.genes <- get.markers(c(6,10))
get.overrepresentation(macrophage.marker.genes)

# SMCs
smc.marker.genes <- get.markers(8)
get.overrepresentation(smc.marker.genes)

# endo
endo.marker.genes <- get.markers(9)
get.overrepresentation(endo.marker.genes)

# b cells
bcells.marker.genes <- get.markers(10)
get.overrepresentation(bcells.marker.genes)

# mast
mast.marker.genes <- get.markers(3)
get.overrepresentation(mast.marker.genes)

# unknown
unknown.marker.genes <- get.markers(4)
get.overrepresentation(unknown.marker.genes)

# plot
fisher.plot <- data.frame(Tcells = c(2.79135,4.953459),
                          Macropgages = c(1.704716, 2.929894),
                          SMCs = c(1.02376, 1.02376),
                          Endothelium = c(1.230978,2.863273),
                          Bcells = c(1.980975, 3.252197),
                          MastCells = c(1.453986, 2.407139),
                          UnknownCells = c(1.550983, 0.08153747),
                          reference = c("Bulk1", "Bulk2"))
                          
fisher.plot <- t(fisher.plot)

fisher.plot <- data.frame(reference = rep(c("Bulk1", "Bulk2"),7),
                          cells = c(rep("Tcells",2), rep("Marcophages",2),rep("SMCs",2), rep("Endothelium",2), rep("Bcells",2), rep("MastCells",2), rep("UnknownCells",2)),
                          value = c(2.79135,4.953459, 1.704716, 2.929894, 1.02376, 1.02376,1.230978,2.863273,1.980975, 3.252197, 1.453986, 2.407139, 1.550983, 0.08153747))

fisher.plot <- data.frame(reference = c("Bulk"),
                          cells = c("Tcells","Marcophages","SMCs", "Endothelium","Bcells","MastCells","UnknownCells"),
                          value = c(2.79135,1.704716, 1.02376, 1.230978,1.980975, 1.453986, 1.550983))                          

ggplot(fisher.plot, aes(fill = reference, x = cells , y = value)) + geom_bar(stat = "identity", position = "dodge") +
  ylab("fold difference") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(paste(result.folder, "/", Sys.Date(), " Fold difference relative expression marker genes pseudobulk vs bulk.pdf", sep = ""))



#------------------------------------------------------------------------------------------------------
# plot all top 10s against Bulk1
# add linear regression line?
plot.pseudo <- log2(single.vs.whole[,c(12,14)] + 0.1)
ggplot(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2])) + 
  geom_point( colour = "black", size = 1) +
  #geom_smooth(method = 'lm', formula=y~x-1, se = F, col = "red") +
  # CD3+CD8+ T-Cells I
  geom_point(data = plot.pseudo[top.10s[,1],], aes(plot.pseudo[top.10s[,1],1], plot.pseudo[top.10s[,1],2]),
             colour = "black", size = 4.2, shape = 25, fill = colour.palette[1]) +
  # CD3+CD4+ T-Cells
  geom_point(data = plot.pseudo[top.10s[,2],], aes(plot.pseudo[top.10s[,2],1], plot.pseudo[top.10s[,2],2]),
             colour = "black", size = 4, shape = 22, fill = colour.palette[2]) +
  # CD3+CD8+ T-Cells III
  geom_point(data = plot.pseudo[top.10s[,6],], aes(plot.pseudo[top.10s[,6],1], plot.pseudo[top.10s[,6],2]),
             colour = "black", size = 3, shape = 24, fill = colour.palette[6] ) +
  # CD3+CD8+ T-Cells II
  geom_point(data = plot.pseudo[top.10s[,4],], aes(plot.pseudo[top.10s[,4],1], plot.pseudo[top.10s[,4],2]),
             colour = "black", size = 3, shape = 23, fill = colour.palette[4]) +
  # ABCA9+LEPR+ Unknown Cells
  geom_point(data = plot.pseudo[top.10s[,3],], aes(plot.pseudo[top.10s[,3],1], plot.pseudo[top.10s[,3],2]),
             colour = "black", shape = 21, fill = colour.palette[3], size = 1.5) +
  # CD14+CD68+ Macrophages I
  geom_point(data = plot.pseudo[top.10s[,5],], aes(plot.pseudo[top.10s[,5],1], plot.pseudo[top.10s[,5],2]),
             colour = "black", shape = 21, fill = colour.palette[5], size = 1.5) +
  # MYH11+ Smooth Muscle Cells
  geom_point(data = plot.pseudo[top.10s[,7],], aes(plot.pseudo[top.10s[,7],1], plot.pseudo[top.10s[,7],2]),
             colour = "black", shape = 21, fill = colour.palette[7], size = 1.5) +
  # CD34+ Endothelial Cells
  geom_point(data = plot.pseudo[top.10s[,8],], aes(plot.pseudo[top.10s[,8],1], plot.pseudo[top.10s[,8],2]),
             colour = "black", shape = 21, fill = colour.palette[8], size = 1.5) +
  # CD14+CD68+ Macrophages II
  geom_point(data = plot.pseudo[top.10s[,9],], aes(plot.pseudo[top.10s[,9],1], plot.pseudo[top.10s[,9],2]),
             colour = "black", shape = 21, fill = colour.palette[9], size = 1.5) +
  # CD79A+ B-Cells
  geom_point(data = plot.pseudo[top.10s[,10],], aes(plot.pseudo[top.10s[,10],1], plot.pseudo[top.10s[,10],2]),
             colour = "black", shape = 21, fill = colour.palette[10], size = 1.5) +
  # KIT+ Mast Cells
  geom_point(data = plot.pseudo[top.10s[,11],], aes(plot.pseudo[top.10s[,11],1], plot.pseudo[top.10s[,11],2]),
             colour = "black", shape = 21, fill = colour.palette[11], size = 1.5) +
  xlab("Bulk") + ylab("Pseudobulk")
 

colour.palette

ggsave(paste(result.folder, "/", Sys.Date(), " All top10 marker genes by colour6.pdf", sep = ""))

ggplot(data = plot.pseudo, aes(plot.pseudo[,1], plot.pseudo[,2])) + 
  geom_point(colour = "black", size = 1) + geom_smooth(method = 'lm')
