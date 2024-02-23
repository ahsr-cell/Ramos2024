## Post genetic demultiplexing - converging scDbl and genetic demultiplexing results 

#Cellranger barcodes
raw_data <- Read10X(data.dir ="/Filers/home/a/andramos/t1data/sc_rnaseq/analysis/cellranger/2feb/D52N_miSFITs/outs/filtered_feature_bc_matrix")
seuobj <- CreateSeuratObject(counts = raw_data, project = "car_miSFITs", min.cells = 3, min.features = 200)

CRBarcodes <- data.frame(colnames(seuobj),stringsAsFactors = FALSE)
names(CRBarcodes)[1] <- "barcode"

#Genetic demultiplexing raw output, see Genetic Demultiplexing
GDBarcodes <- read.delim("/project/tsslab/shared/ramoscandido/6_tenxpreprocessing/10X-demux/gendem/gPlexA/results.gPlexA.dir/all_demultiplexing_results.tsv")
names(GDBarcodes)[1] <- "barcode"
names(GDBarcodes)[2] <- "aware"
names(GDBarcodes)[3] <- "blind"
GDBarcodes$barcode <- gsub("-gPlexA","",as.character(GDBarcodes$barcode))
GDBarcodes$aware <- gsub("UNASSIGNED","DOUBLET",as.character(GDBarcodes$aware)) #Filter out unassigned, mark as doublet
GDBarcodes$blind <- gsub("UNASSIGNED","DOUBLET",as.character(GDBarcodes$blind)) #filter out unassigned, mark as doublet
GDBarcodes$aware <- gsub("18-AD2-4-N704", "",as.character(GDBarcodes$aware))
GDBarcodes$aware <- gsub("18-AD2-1-N701", "",as.character(GDBarcodes$aware))
GDBarcodes$aware <- gsub("18-AD2-2-N702", "",as.character(GDBarcodes$aware))

GDBarcodes_filter_aware <- filter( GDBarcodes, aware == "DOUBLET") 
GDBarcodes_filter_blind <- filter( GDBarcodes, blind == "DOUBLET") 
GDBarcodes_filter_out <- full_join(GDBarcodes_filter_aware, GDBarcodes_filter_blind, "barcode")
GDBarcodes_filter_out$ids <- "DOUBLET"
GDBarcodes_filter_out <- GDBarcodes_filter_out[c(1,6)]
GDBarcodes_keep <- anti_join(GDBarcodes, GDBarcodes_filter_out, "barcode")
GDBarcodes_keep$ids <- GDBarcodes_keep$aware
GDBarcodes_keep <- GDBarcodes_keep[c(1,4)]

GDBarcodesFinal <- rbind(GDBarcodes_keep, GDBarcodes_filter_out)

#scDblFinder barcodes: in silico check if genetic demultiplexing missed any doublets - generated in 0_Step0_EmptyDroplets_DoubletsRemoval.R, line 51
names(scDblBarcodes)[1] <- "barcode"
names(scDblBarcodes)[2] <- "ids"
scDblDoublets <- filter(scDblBarcodes, ids == "doublet")  
CheckDbls <- anti_join(scDblDoublets, GDBarcodesFinal, "barcode") #results in an empty df, therefore GD covered all doublets

Barcodes <- semi_join(GDBarcodesFinal, CRBarcodes, "barcode") #merge 

write.table(Barcodes, file = "GDscDbl_annotated_Barcodes.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE) ## used in next step

