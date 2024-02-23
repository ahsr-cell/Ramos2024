###SCENIC: Genie3 on Cluster - hg38 - with both 100/500bp and 10kb

require(Seurat)
require(GENIE3)
require(SingleCellExperiment)
require(AUCell)
require(RcisTarget)
require(SCopeLoomR)
require(KernSmooth)
require(BiocParallel)
require(ggplot2)
require(data.table)
require(grid)
require(ComplexHeatmap)

knitr::opts_knit$set(root.dir="home/a/andramos/t1data/sc_rnaseq/analysis/R_analysis/Analysis/SCENIC/hg38/refseq") 
print(paste(Sys.time(), "set working directory"))

load("/home/a/andramos/t1data/sc_rnaseq/analysis/R_analysis/Analysis/SCENIC/hg38/PreProcessing_for_hg38SCENIC_18July.RData")
exprMat <- FetchData(subsetPR, vars = rownames(subsetPR))
exprMat <- as.data.frame(t(exprMat), stringsAsFactors = FALSE)
cellInfo <- data.frame(subsetPR$Percentile_Rank, stringsAsFactors = FALSE)
names(cellInfo)[1] <- "Percentile_Rank"
exprMat <- as.matrix(exprMat)
print(paste(Sys.time(), "loaded in prepared subsetPR, exprMat, cellInfo"))

library(SCENIC)
org <- "hgnc" 
dbDir <- "/Users/ahsramos/Desktop/Dissertation/work_in_progress/scRNAseq_bioinformatic_analysis/SCENIC/hg38/pySCENIC/cisTargetDatabase" # RcisTarget databases location
myDatasetTitle <- "D52NCAR_PR_500_100_10kb" # choose a name for your analysis
dbs <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
scenicOptions <- initializeScenic(org="hgnc", dbDir = "/home/a/andramos/t1data/sc_rnaseq/analysis/R_analysis/Analysis/SCENIC/cisTarget_databases/hg38" , dbs=c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"), datasetTitle=myDatasetTitle, nCores=8) 
scenicOptions@settings$db_mcVersion <- "v9"
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
print(paste(Sys.time(), "SCENIC options created"))

genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

print(paste(Sys.time(), "running correlation"))
runCorrelation(exprMat_filtered, scenicOptions)
print(paste(Sys.time(), "correlation completed, starting genie3 run"))

#exprMat_filtered_log <- log2(exprMat_filtered+1) - I inputted a log-normalised, filtered seuobj 
runGenie3(exprMat_filtered, scenicOptions)
print(paste(Sys.time(), "success-genie3 completed"))
