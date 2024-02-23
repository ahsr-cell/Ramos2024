#A script for demultiplexing 10X cell barcodes and miSFITs post dataframe creation
library(tidyverse)
library(ggbreak) 

#import dataframe of counts
## provided - see 0_Step2_misfit_demultiplexing_all10xcbs.R for creation
misdm <- read.delim("/Users/ahsramos/Desktop/Dissertation/work_in_progress/scRNAseq_bioinformatic_analysis/Demultiplexing_Statistics/misfits_demultiplexing_all10xcbs_step1.txt")

#for indexing purposes, grab the barcodes column
tenxcellbs <- misdm$barcodes
#as max() functions do not work with characters, subset misdm_df to just read counts, use this for evaluating counts
working_misdm <- dplyr::select(misdm,-c(1,))
#define minimum read count classification
min_rc = 25
#define minimum factor classification
min_fc = 3
#create new column to store results
misdm$Results <- NA
misdm$HighestCount <- NA
misdm$SecondHighestCount <- NA
misdm$Classification <- NA 

#starting loop: go through cell barcode through cell barcode
for (i in 1:length(tenxcellbs)) {
  #pull current cell barcode counts
  z <- working_misdm %>% dplyr::slice(i)
  #within that row, identify the highest count
  x <- max(z, na.rm = TRUE)
  misdm[i,23] <- x
  #within that row, identify the second highest count
  y <- max(z[z != max(z, na.rm = TRUE)], na.rm = TRUE) 
  misdm[i,24] <- y
  #progress report
  print(paste(Sys.time(), "progress is",i, "out of", length(tenxcellbs)))
  #classification 1: does the highest miSFIT read count equal the second highest read count (e.g. 1 vs 1, 0 vs 0)
  if ( x == y) {
    misdm[i,22] <- "NA"
    misdm[i,25] <- "Highest equals second highest"
    print("10X cell barcode assigned NA - highest miSFIT read count equals second highest misFIT read count")
    #classification 2: does the highest miSFIT read count meet the minimum read count
  } else if ( x <= min_rc) {
    misdm[i,22] <- "NA"
    misdm[i,25] <- "Below mrc"
    print("10X cell barcode assigned NA - highest miSFIT read count below minimum read count threshold")
    #classification 3: does the highest miSFIT read count meet the minimum factor?
  } else if ( x <= y*min_fc) {
    misdm[i,22] <- "NA"
    misdm[i,25] <- "Fails mfc"
    print("10X cell barcode assigned NA - highest miSFIT read count does not meet minimum factor thresold")
  } else {
    max_pos<- matrix(c(1, apply(z, 1, which.max)), ncol=2)
    #store column index
    max_pos[,2] <- max_pos[,2] + 1
    #+1 to account for removal of barcodes column
    misdm[i,22] <- colnames(misdm)[max_pos[,2]]
    misdm[i,25] <- "Assigned"
    #in column of current row index, write miSFIT barcode, using stored index
    print("miSFIT barcode assigned to 10X cell barcode")
  }
}
#max() function will interpret 0 as second highest as -Inf, therefore let's replace it with 0
misdm[misdm == -Inf] <- 0

#Summary stats 
miSFIT_dm <- misdm[,c("barcodes","Results")] ### used in 1_Startup - saved file as misfit_demultiplexing_10xbarcodes_complete.txt
miSFIT_dm$barcodes <- paste0(miSFIT_dm$barcodes, "-1")
summary_stats <- misdm[,c("barcodes", "Results", "HighestCount","SecondHighestCount")]

###Filtering demultiplexing from 18,761 barcodes to 13,921 barcodes 
load("~/Desktop/Dissertation/work_in_progress/scRNAseq_bioinformatic_analysis/scRNAseq_baseforlocalanalysis.RData")
#creating seuratobj
seuobj <- CreateSeuratObject(counts = raw_data, project = "car_miSFITs", min.cells = 3, min.features = 200)
#match barcodes info to cellranger barcodes and add to metadata
myBarcode = rownames(seuobj@meta.data)
ids = barcodes[match(myBarcode,barcodes$barcode),]
seuobj$ids = ids$id
#adding in mitochondrial percentage to metadata
seuobj[["percent.mt"]] <-PercentageFeatureSet(seuobj, pattern = "^MT-")
#filter and subset dataset by that 10%> percent.mt AND has had its donor identified 
subset <- subset(seuobj, 
                 subset = nFeature_RNA >200 & 
                   nFeature_RNA <6000 & percent.mt <10 & 
                   (ids=="A" | ids=="B" | ids=="C"))
bc13921 <- colnames(subset)
bc13921 <- as.data.frame(bc13921, stringsAsFactors = FALSE)
colnames(bc13921)[1] <- "barcodes"
bc13921$barcodes <- gsub("-1","",as.character(bc13921$barcodes))
misdm_13921 <- semi_join(misdm, bc13921, by = "barcodes")
summary_stats <- misdm_13921[,c("barcodes", "Results", "HighestCount","SecondHighestCount")]

#Extended Data Figure 12b
ggplot(summary_stats, aes(x=HighestCount)) +
  theme_minimal() +
  geom_histogram(binwidth=20, fill = "steelblue") +
  scale_y_break(c(1000, 8500)) + 
  xlab("Highest count") +
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
ggplot(summary_stats, aes(x=HighestCount)) +
  theme_minimal() +
  geom_histogram(binwidth=5, fill = "steelblue") +
  xlab("Highest count") +
  xlim(1,1000) +
  geom_vline(xintercept = 25, linetype="dashed", color = "red") +
  annotate("text", x = 200, y = 162.5, label = "minimum read count = 25", color = "red") +  
  ylim(0,950) +
  scale_y_break(c(200, 900)) + 
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

#Extended Data Figure 12c
summary_stats$Ratio <- mapply('/', summary_stats$HighestCount, summary_stats$SecondHighestCount)
Ratio <- summary_stats %>% dplyr::filter(Ratio != "NaN", Ratio != "Inf")
ggplot(summary_stats, aes(x=Ratio)) +
  geom_histogram(binwidth=0.5, fill = "steelblue") + 
  theme_minimal() + 
  xlab("Ratio") +
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggplot(Ratio, aes(x=Ratio)) +
  geom_histogram(binwidth=0.3, fill = "steelblue") + 
  xlab("Ratio") +
  xlim(0,50) +
  ylim(0,500) +
  theme_minimal() + 
  geom_vline(xintercept = 3, linetype="dashed", color = "red") +
  annotate("text", x = 9.5, y = 425, label = "minimum factor = 3", color = "red") +  
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

##Number of called vs uncalled
a <- as.data.frame(c("Called", "Uncalled"))
a$Values <- c((length(Called$Results)),(length(misdm_13921$Results) - length(Called$Results)))
names(a)[1] <- "Classification"

#Extended Data Figure 12d
ggplot(a, aes(x=Classification, y=Values)) + 
  geom_bar(stat="identity",fill="steelblue", color ="black", size = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  geom_text(aes(label=Values), size=5, vjust = 1.5, colour = "white") +
  ylab("Number") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text=element_text(size=10),
        axis.title.x = element_blank())

###Barplot of #cells/miSFIT
cells_per_misfit <- as.data.frame(c("Perfect1x", "Perfect2x","Perfect4x","Scramble","Empty",
             "V1","V2","V3","V4","V5",
             "V6", "V7","V8","V9","V10",
             "V11","V12","V13","V14","V15"))
#, "NA"
cells_per_misfit$Number_of_cells <- c(sum(summary_stats$Results == "CCTACCTGCACTGTAAGCACTTTGA"),           
                  sum(summary_stats$Results == "CTACCTGCACTGTAAGCACTTTGTAT"),           
                  sum(summary_stats$Results == "CTACCTGCACTGTAAGCACTTTGTATCTACCTGCACTGTAAGCACTTTGTAT"),           
                  sum(summary_stats$Results == "CGCGCTTCCGCGGCCCGTTCAAG"),           
                  sum(summary_stats$Results == "CGTAGGCGCGCCGTCTCTACG"),           
                  sum(summary_stats$Results == "CTACCTGAACTGTAAGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCAGTGTACGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCTCTGTAAGCACTTTG"),           
                  sum(summary_stats$Results == "CAACCTGCACTGTAAGTACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCGCTGTAAGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCACTGTAAGCACTTGG"),           
                  sum(summary_stats$Results == "ATACATGCACTGTAAGCACTTTG"),           
                  sum(summary_stats$Results == "CTTCCTGCACTGAAAGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCACTGTAAGCGCTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCACTCTAAGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCGCTGTACGCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCACTCTAAGCACTGTG"),           
                  sum(summary_stats$Results == "CTACCCGCACTGTAACCACTTTG"),           
                  sum(summary_stats$Results == "CTACCTGCACTGTAAGAGCTTTG"),           
                  sum(summary_stats$Results == "CTCCCTGCACTGTAAGCACTTTG"))
#sum(summary_stats$Results == "NA")
names(cells_per_misfit)[1] = "miSFIT"
level_order <- c("Perfect1x", "Perfect2x","Perfect4x","Scramble","Empty",
                 "V1","V2","V3","V4","V5",
                 "V6", "V7","V8","V9","V10",
                 "V11","V12","V13","V14","V15")

cells_per_misfit <- cells_per_misfit[order( cells_per_misfit[,2] ),]
cellspm_order <- cells_per_misfit[,1]
#,"NA"

#Extended Data Figure 12e
ggplot(cells_per_misfit, aes(x=factor(miSFIT, level = cellspm_order), y=Number_of_cells)) + 
  geom_bar(stat="identity",colour = "black", fill="steelblue", size = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  geom_text(aes(label=Number_of_cells), size=3, vjust = -0.25, colour = "black") +
  xlab("miSFIT") +
  ylab("Number of cells") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

