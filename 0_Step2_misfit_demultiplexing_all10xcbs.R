## R Script for miSFIT demultiplexing
rm(list=ls(all=TRUE))
setwd("/home/a/andramos/t1data/sc_rnaseq/analysis/R_analysis/miSFIT_matching")

library(ShortRead)
library(tidyverse)
library(Seurat)
library(reshape)
library(Biostrings)
print(paste(Sys.time(), "Finished loading packages, loading seuratobj creation"))

###Loading in code to attain list of barcodes that incorporate cellranger QC filtering, Seurat QC filtering, and gendem+doubletremoval 
#cellranger filtering
raw_data <- Read10X(data.dir ="/Filers/home/a/andramos/t1data/sc_rnaseq/analysis/cellranger/2feb/D52N_miSFITs/outs/filtered_feature_bc_matrix")
#seuratobj creation
seuobj <- CreateSeuratObject(counts = raw_data, project = "car_miSFITs", min.cells = 3, min.features = 200)
#add in mitochondrial percentage to metadata
seuobj[["percent.mt"]] <-PercentageFeatureSet(seuobj, pattern = "^MT-")
print(paste(Sys.time(), "Finished seuratobj creation, importing gendem barcodes"))

#loading in genetic demultiplexing, doublet, and empty droplet filtering results 
## provided - see 0_Step1_DemultiplexingConverge.R for creation
#note: if to be filtered out, labeled as DOUBLET
barcodes <- read.delim("/t1-data/project/tsslab/andramos/sc_rnaseq/analysis/R_analysis/GDscDbl_annotated_Barcodes.txt")

#match barcodes info to cellranger barcodes and add to metadata
myBarcode = rownames(seuobj@meta.data)
ids = barcodes[match(myBarcode,barcodes$barcode),]
seuobj$ids = ids$id
print(paste(Sys.time(), "Added gendem barcodes to seuobj, creating subset"))

#filter and subset dataset by that 10%> percent.mt 
subset <- subset(seuobj, 
                 subset = nFeature_RNA >200 & 
                   nFeature_RNA <6000 & percent.mt <10)  
#                   & (ids=="A" | ids=="B" | ids=="C")) 
print(paste(Sys.time(), "Finished subset, extracting 10x cell barcodes for misfit demultiplexing"))

#extracting post QC barcodes and creating a dataframe to serve as a basis for recording
wbcs <- subset@active.ident
wbcs <- as.data.frame(wbcs)
wbcs <- tibble::rownames_to_column(wbcs, "barcodes")
wbcs <- wbcs[-c(2)]
wbcs$barcodes<-gsub("-1","",as.character(wbcs$barcodes))
print(paste(Sys.time(), "Creating template dataframe to serve as recording basis, adding misfit barcodes"))

#wbcs <- magrittr::set_rownames(wbcs, wbcs$barcodes)
#adding in each miSFIT barcode
#perfect1x
wbcs$CCTACCTGCACTGTAAGCACTTTGA <- NA
#perfect2x
wbcs$CTACCTGCACTGTAAGCACTTTGTAT <- NA
#perfect4x
wbcs$CTACCTGCACTGTAAGCACTTTGTATCTACCTGCACTGTAAGCACTTTGTAT <- NA
#scr
wbcs$CGCGCTTCCGCGGCCCGTTCAAG <- NA
#empty
wbcs$CGTAGGCGCGCCGTCTCTACG <- NA
#V1
wbcs$CTACCTGAACTGTAAGCACTTTG <- NA
#V2
wbcs$CTACCTGCAGTGTACGCACTTTG <- NA
#V3
wbcs$CTACCTGCTCTGTAAGCACTTTG <- NA
#V4
wbcs$CAACCTGCACTGTAAGTACTTTG <- NA
#V5
wbcs$CTACCTGCGCTGTAAGCACTTTG <- NA
#V6
wbcs$CTACCTGCACTGTAAGCACTTGG <- NA
#V7
wbcs$ATACATGCACTGTAAGCACTTTG <- NA
#V8
wbcs$CTTCCTGCACTGAAAGCACTTTG <- NA
#V9
wbcs$CTACCTGCACTGTAAGCGCTTTG <- NA
#V10
wbcs$CTACCTGCACTCTAAGCACTTTG <- NA
#V11
wbcs$CTACCTGCGCTGTACGCACTTTG <- NA
#V12
wbcs$CTACCTGCACTCTAAGCACTGTG <- NA
#V13
wbcs$CTACCCGCACTGTAACCACTTTG <- NA
#V14
wbcs$CTACCTGCACTGTAAGAGCTTTG <- NA
#V15
wbcs$CTCCCTGCACTGTAAGCACTTTG <- NA
#wbcs <- wbcs[-c(1)]
print(paste(Sys.time(), "Finished adding misfit barcodes as columns, beginning to grab misfit demultiplexing fastq files, starting with read1"))

#Grabbing R1 read for 10Xcellbarcode
r1<- readDNAStringSet("/home/a/andramos/t1data/sc_rnaseq/data/10x_misfits_raw_fastq/fastq_gz/R1_cellbc/i7andi5-Indexing-PairedEnd_S1_L001_R1_001.fastq",format="fastq")
r1<-as.data.frame(r1)
r1 <- tibble::rownames_to_column(r1, "identifier")
#r1$identifier<-gsub("M01913:28:000000000-JVGK9:1:1101:","",as.character(r1$identifier))
#r1$identifier<-gsub(" 1:N:0:1","",as.character(r1$identifier))
names(r1)[2] <- "r1"
print(paste(Sys.time(), "Finished R1, 10X cell barcode fastq file import and dataframe creation, importing read2"))

#Grabbing R2 read for misFIT barcode
r2<- readDNAStringSet("/home/a/andramos/t1data/sc_rnaseq/data/10x_misfits_raw_fastq/fastq_gz/R2_misfbc/i7andi5-Indexing-PairedEnd_S1_L001_R2_001.fastq",format="fastq")
r2 <- as.data.frame(r2)
r2 <- tibble::rownames_to_column(r2, "identifier")
#r2$identifier<-gsub("M01913:28:000000000-JVGK9:1:1101:","",as.character(r2$identifier))
#r2$identifier<-gsub(" 2:N:0:1","",as.character(r2$identifier))
names(r2)[2] <- "r2"
print(paste(Sys.time(), "Finished R2, misfit barcode fastq file import and dataframe creation, beginning to merge the two dataframes"))

#creating reads dataframe
reads <- dplyr::inner_join(r1, r2, by = "identifier")
print(paste(Sys.time(), "Creating lists to serve for searching function"))

#extracting barcodes from wbcs dataframe for searching
cell_barcodes <- wbcs$barcodes
misfit <- tail(colnames(wbcs),-1)

print(paste(Sys.time(), "Starting for loops to begin searching and counting"))
# Iterate over cell barcodes
for (i in seq_along(cell_barcodes)) {
  
  # Loop over R1 file till the line is the empty vector.
  print(paste(Sys.time(), "progress is" ,i, "out of",length(cell_barcodes)))
  subset <- dplyr::filter(reads, grepl(cell_barcodes[[i]], r1, ignore.case = TRUE))
  
  for (j in seq_along(misfit)) {
    #count how many instances 
    count <- nrow(dplyr::filter(subset, grepl(misfit[[j]],r2)))
    #tell it where to store it within wbcs
    col <- misfit[[j]]
    wbcs[,col][wbcs$barcodes == cell_barcodes[[i]]] <- count
  }
  
}
print(paste(Sys.time(), "Finished search and count, creating table of output"))
write.table(wbcs,"misfits_demultiplexing_all10xcbs_step1.txt",sep="\t") ###used in 0_Step3
print(paste(Sys.time(), "Script successfully finished"))
