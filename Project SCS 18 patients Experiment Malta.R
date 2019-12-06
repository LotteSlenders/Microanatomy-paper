### 2019-5-16

# Project SCS 18 patients Experiment Malta

# erdmann  paper with single cell data

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("C:/Users/lslender/Documents/Projecten/Project SCS 18 patients/Experiment Malta")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 3.5.3

### load packages
library(Seurat)         # V2.3.4
library(plyr)           # v0.7.8
library(SingleR)        # v0.2.2
library(org.Hs.eg.db)

# object
load(file = "all.seur.combined.rdata")
seuset <- all.seur.combined
rm(all.seur.combined)

# set seed
set.seed(435)

#-----------------------------------------------------------------------------------------
# lipid metabolism
# lipid.metabolism <- c("APOA5", "APOC3", "LPL", "ANGPTL4", "APOC", "LDLR", "PCSK9", "APOB", "ABCG5", "ABCG8", "HNF1A", "ABO",
#                       "TRIB1", "SORT1", "LPA", "LRP1", "SCARB1", "PLTP")
lipid.metabolism <- c("ANGPTL4", "LDLR", "PCSK9", "APOB", "ABCG5", "ABCG8", "HNF1A", "ABO",
                      "TRIB1", "SORT1", "LPA", "LRP1", "SCARB1", "PLTP")

lapply(lipid.metabolism, function(x){
  VlnPlot(seuset, features.plot = x)
})

# blood pressure
blood.pressure <- c("NOS3", "SH2B3", "CYP17A1", "GUY1A1", "FURIN", "AGT", "ARHGAP42")

# Mitosis proliferation
mitosis.proliferation <- c("ZC3HC1", "CDKN2A", "RHOA", "MAD2L1", "MAD1L1", "CFDP1", "BCAS3", "MRAS", "RAD50",
                           "STAG1", "CENPW", "CDC123", "PDS5B", "COPRS", "SMARCA4", "RHOA")

# nevascularization angiogenesis
neovascularization <- c("DAB2IP", "SMAD3", "FGD6", "ANKS1A", "BCAS3", "VEGFA", "TGFB1", "CCM2", "ZFPM2")

# NO signalling
no.signalling <- c("GUCY1A1", "EDNRA", "NOS3", "EDN1", "PDE5A", "PDE3A", "TBAS1", "ARHGAP42")

# Vascular remodeling
vascular.remodeling <- c("COL4A1", "COL4A2", "COL6A3", "MIA3", "REST-NOA1", "TSPAN14", "PDGFD", "SWAP70", "KSR2",
                         "ADAMTS7", "BCAS3", "FLT1", "PECAM1", "SH3PXD2A", "PI16", "SERPINH1", "LOX", "ITGB5",
                         "HTRA1", "FN1", "TSN1", "FURIN", "BMP1", "RPL17")

# transcriptome gene regulation
transcriptome.regulation <- c("CTR9", "PMAIP1", "TDRKH", "PRIM2", "FHL3", "YY1", "FOXC1", "MAP3K1", "KLF4", "DAB2IP",
                              "ZNF589", "RGS12", "HDAC9", "ARNTL", "FGF5", "N4BP2L2", "BACH1", "ZNF831", "PCIF1",
                              "HNRNPD", "SKI")

# Inflammation
inflammation <- c("IL5", "CXCL12", "MRAS", "PLG", "TRIM22", "C1S", "SH3B3", "IL6R", "FAM213A", "PRKCE", "CFTR", "DHX58",
                  "TRIM5", "TRIM6")

# unknown 
# no HDGFL1
unknown1 <- c("ATP1B1", "NME7", "DDX59", "CAMSAP2", "LMOD1", "TEX41", "NBEAL1", "IRS1", "PLCG1", "ZNF827", "SLC22A4", "ARHGAP26",
              "BCAP29", "GPR22", "KLHDC10", "PARP12", "FNDC3B")
# no STDB1, HP
unknown2 <- c("DNAJC13", "NCK1", "PPP2R3A", "PALLD", "TIPARP", "FIGN", "HHAT", "PRDM16", "HHIPL1", "TRIP4", "MFGE8",
              "CFDP1", "BCAR1", "CDH13", "SMG6", "RNASE13")
# no KCNE2
unknown3 <- c("ZNF507", "HNRNPUL1", "SNRPD2", "PROCR", "EIF6", "ADORA2A", "SHROOM3", "SERPINA1", "PLCG1", "FCHO1",
              "HSD17B12", "PSMA3", "MCF2L", "SIPA1", "PRDM8")
unknown4 <- c("CORO6", "ANKRD13B", "NDUFA12", "MAP3K7CL", "RAC1", "NAT2", "HOXC4", "CCDC92", "CDKN1A", "PRIM2", "PLEKHG1",
              "PNMT", "GOSR2", "UBE2Z", "DDX5", "CCDC92")
unknown5 <- c("UNC5C", "ALS2CL", "RTP3", "GIP", "TEX2", "MRVI1", "ARHGEF26", "KCNJ13")


unknown <- c(unknown1, unknown2, unknown3, unknown4, unknown5)
length(unknown)
unknown <- unique(unknown)
length(unknown)

#-----------------------------------------------------------------------------------------
### prep data
## get average seuset
average.seuset <- AverageExpression(seuset, return.seurat = T)

## do k means
seed <- 63

Kmeans.object <- DoKMeans(average.seuset,
                          genes.use = unknown,
                          k.genes = 3,
                          do.plot = T,
                          k.seed = seed)
ggsave(paste(result.folder, "/", Sys.Date(), " Unkown genes Erdmann paper vs scRNAseq2.pdf", sep = ""))

#-----------------------------------------------------------------------------------------
DEG <- read.delim("All_DEG_distributed_across_Kmeans_12062019.txt", stringsAsFactors = F, sep = " ")

unknown %in% rownames(DEG)
g <- unknown[unknown %in% rownames(DEG)]
DEG$ident <- rownames(DEG)
View(DEG[g,])


#-----------------------------------------------------------------------------------------
