#A function to calculate the frequency of each miSFIT17 variant in the transduction pool

#cleaning up
rm(list=ls(all=TRUE)) 

#get set up: work here, find my files and functions
setwd("~/Desktop/MiSeq_2/3_Raw")
in_file_list <- list.files()
source("~/Desktop/R_Scripts/miSFITgDNA_cDNA_Counts_Functions.R")

#Load packages
library(Biostrings)
library(stringr)

#These are my variants, compile all of them into a list
ref_perfect = DNAString("ggTcaaagtgcttacagtgcaggtagGCT")
ref_perfect2x <- DNAString("ggTcaaagtgcttacagtgcaggtagATAcaaagtgcttacagtgcaggtagGCT")
ref_perfect4x <- DNAString("caaagtgcttacagtgcaggtagATAcaaagtgcttacagtgcaggtagATA")
ref_scramble_1 <- DNAString("CTTGAACGGGCCGCGGAAGCGCG")
ref_empty <- DNAString("TCGTAGAGACGGCGCGCCTACG")
ref_v1 <- DNAString("caaagtgcttacagtTcaggtag")
ref_v2 <- DNAString("caaagtgcGtacaCtgcaggtag")
ref_v3 <- DNAString("caaagtgcttacagAgcaggtag")
ref_v4 <- DNAString("caaagtActtacagtgcaggtTg")
ref_v5 <- DNAString("caaagtgcttacagCgcaggtag")
ref_v6 <- DNAString("cCaagtgcttacagtgcaggtag")
ref_v7 <- DNAString("caaagtgcttacagtgcaTgtaT")
ref_v8 <- DNAString("caaagtgcttTcagtgcaggAag")
ref_v9 <- DNAString("caaagCgcttacagtgcaggtag")
ref_v10 <- DNAString("caaagtgcttaGagtgcaggtag")
ref_v11 <- DNAString("caaagtgcGtacagCgcaggtag")
ref_v12 <- DNAString("caCagtgcttaGagtgcaggtag")
ref_v13 <- DNAString("caaagtgGttacagtgcGggtag")
ref_v14 <- DNAString("caaagCTcttacagtgcaggtag")
ref_v15 <- DNAString("caaagtgcttacagtgcaggGag")
miSFITs <- list(ref_perfect, ref_perfect2x, ref_perfect4x, ref_scramble_1, ref_empty, ref_v1, ref_v2, ref_v3, ref_v4, ref_v5, ref_v6, ref_v7, ref_v8, ref_v9, ref_v10, ref_v11, ref_v12, ref_v13, ref_v14, ref_v15)

#In my workspace are my fastq files, on an indice/individual-basis import them as DNAStrings
#For the current fastq file imported, utilise my function to calculate the frequencies and absolute counts for each variant listed in my compiled variants list
#Put the output of these functions into respectively labeled holders
#Take these holders and compile them into a single place
#Save this data in the current working directory
for(i in 1:length(in_file_list))
{
  print(paste("current file is ",in_file_list[i]))
cgDNAs<- readDNAStringSet(in_file_list[i],format="fastq")
total_reads<- length(cgDNAs)
Freqs_Holder<-Frequencies(miSFITs,cgDNAs)
Absols_Holder<-Absolutes(miSFITs,cgDNAs)
results_list<-list(total_reads,miSFITs,Freqs_Holder,Absols_Holder)
save(results_list,file=paste(in_file_list[i],"_validation.Rdata",sep=""))
}

