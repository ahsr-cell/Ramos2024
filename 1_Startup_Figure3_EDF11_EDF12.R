require(Seurat)
require(tidyverse)
require(reshape)
require(viridis)
require(ggbreak)
require(ggrepel)
require(drc)
require(SingleCellExperiment)
require(infotheo)
require(tidyinftheo)
require(svglite)
require(pheatmap)

###Creation of Seurat Object, data wrangling, and annotating cells by miSFIT and CAR expression percentile rank
#loading in cellranger results
raw_data <- Read10X(data.dir ="//Filers/home/a/andramos/t1data/sc_rnaseq/analysis/cellranger/2feb/D52N_miSFITs/outs/filtered_feature_bc_matrix")
#creating seuratobj
seuobj <- CreateSeuratObject(counts = raw_data, project = "car_miSFITs", min.cells = 3, min.features = 200)

#loading in genetic demultiplexing, doublet, and empty droplet filtering results
## provided - see 0_Step1_DemultiplexingConverge.R for creation
#note: if to be filtered out, labeled as DOUBLET
barcodes <- read.delim("/t1-data/project/tsslab/andramos/sc_rnaseq/analysis/R_analysis/GDscDbl_annotated_Barcodes.txt")

#loading in misfit demultiplexing results, provided - see 0_Step3_miSFIT_demultiplexing_processing.R for creation
miSFIT_dm <- read.csv("/t1-data/project/tsslab/andramos/sc_rnaseq/analysis/R_analysis/Analysis/miSFIT_matching/misfit_demultiplexing_10xbarcodes_complete.txt")

#For genetic demultiplexing: match barcodes info to cellranger barcodes and add to metadata
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

#For miSFIT demultiplexing: match misfits and barcodes info to cellranger barcodes and add to metadata
myBarcode_a = rownames(subset@meta.data)
miSFITs = miSFIT_dm[match(myBarcode_a,miSFIT_dm$barcodes),]
subset$miSFITs = miSFITs$Results
#normalise
subset <- NormalizeData(subset, normalization.method = "LogNormalize", scale.factor = 10000)
#id variable features
subset <- FindVariableFeatures(subset, selection.method = "vst")
#view variable features
top50sub <- head(VariableFeatures(subset), 50)
#scale data
all.genes.subset <- rownames(subset)
subset <- ScaleData(subset, features=all.genes.subset)
#PCA
subset <- RunPCA(subset, features = VariableFeatures(object = subset))
#adding in cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
subset <- CellCycleScoring(subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
subset$CC.Difference <- subset$S.Score - subset$G2M.Score

#Neighbor and cluster identification
subset <- FindNeighbors(subset, dims=1:50)
subset <- FindClusters(subset, resolution = 0.5)
#UMAP & TSNE
subset <- RunUMAP(subset, dims = 1:50)
subset <- RunTSNE(subset, dims = 1:50)
#Adding to metadata
subset[["Clustering"]] <-Idents(object = subset)
Idents(object = subset) <- "Clustering"
subset@active.ident <- factor(subset@active.ident, levels=c("4","1","8","3","7","6","5","0","2")) 
subset$Clustering <- factor(subset$Clustering,
                            levels = c("4", "1", "8", "3", "7", "6", "5", "0", "2"))

cluster0 <- WhichCells(subset, idents=0)
cluster1 <- WhichCells(subset, idents=1)
cluster2 <- WhichCells(subset, idents=2)
cluster3 <- WhichCells(subset, idents=3)
cluster4 <- WhichCells(subset, idents=4)
cluster5 <- WhichCells(subset, idents=5)
cluster6 <- WhichCells(subset, idents=6)
cluster7 <- WhichCells(subset, idents=7)
cluster8 <- WhichCells(subset, idents=8)

#converting to dataframes
cluster0 <- as.data.frame(cluster0)
cluster1 <- as.data.frame(cluster1)
cluster2 <- as.data.frame(cluster2)
cluster3 <- as.data.frame(cluster3)
cluster4 <- as.data.frame(cluster4)
cluster5 <- as.data.frame(cluster5)
cluster6 <- as.data.frame(cluster6)
cluster7 <- as.data.frame(cluster7)
cluster8 <- as.data.frame(cluster8)

#changing column name to use dataframe as basis for filtering
names(cluster0)[1] <- "barcode"
names(cluster1)[1] <- "barcode"
names(cluster2)[1] <- "barcode"
names(cluster3)[1] <- "barcode"
names(cluster4)[1] <- "barcode"
names(cluster5)[1] <- "barcode"
names(cluster6)[1] <- "barcode"
names(cluster7)[1] <- "barcode"
names(cluster8)[1] <- "barcode"

#for plotting frequency, let's add a second column with a label for each cluster
cluster0$cluster <- "0"
cluster1$cluster <- "1"
cluster2$cluster <- "2"
cluster3$cluster <- "3"
cluster4$cluster <- "4"
cluster5$cluster <- "5"
cluster6$cluster <- "6"
cluster7$cluster <- "7"
cluster8$cluster <- "8"

#let's conglomerate these 'sub'-dfs into one single one
clusters <- dplyr::bind_rows(cluster0, cluster1)
clusters <- dplyr::bind_rows(clusters, cluster2)
clusters <- dplyr::bind_rows(clusters, cluster3)
clusters <- dplyr::bind_rows(clusters, cluster4)
clusters <- dplyr::bind_rows(clusters, cluster5)
clusters <- dplyr::bind_rows(clusters, cluster6)
clusters <- dplyr::bind_rows(clusters, cluster7)
clusters <- dplyr::bind_rows(clusters, cluster8)

#Binning by miSFITs
Idents(object = subset) <- "miSFITs"
perf1x <- as.data.frame(WhichCells(subset, idents = c("CCTACCTGCACTGTAAGCACTTTGA")))
perf2x <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTGTAAGCACTTTGTAT")))
#perf4x <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTGTAAGCACTTTGTATCTACCTGCACTGTAAGCACTTTGTAT")))
scr <- as.data.frame(WhichCells(subset, idents = c("CGCGCTTCCGCGGCCCGTTCAAG")))
empty <- as.data.frame(WhichCells(subset, idents = c("CGTAGGCGCGCCGTCTCTACG")))
v1 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGAACTGTAAGCACTTTG")))
v2 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCAGTGTACGCACTTTG")))
v3 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCTCTGTAAGCACTTTG")))   
v4 <- as.data.frame(WhichCells(subset, idents = c("CAACCTGCACTGTAAGTACTTTG")))  
v5 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCGCTGTAAGCACTTTG")))  
v6 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTGTAAGCACTTGG")))
v7 <- as.data.frame(WhichCells(subset, idents = c("ATACATGCACTGTAAGCACTTTG"))) 
v8 <- as.data.frame(WhichCells(subset, idents = c("CTTCCTGCACTGAAAGCACTTTG"))) 
v9 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTGTAAGCGCTTTG")))
v10 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTCTAAGCACTTTG")))
v11 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCGCTGTACGCACTTTG")))
v12 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTCTAAGCACTGTG")))
v13 <- as.data.frame(WhichCells(subset, idents = c("CTACCCGCACTGTAACCACTTTG")))
v14 <- as.data.frame(WhichCells(subset, idents = c("CTACCTGCACTGTAAGAGCTTTG")))  
v15 <- as.data.frame(WhichCells(subset, idents = c("CTCCCTGCACTGTAAGCACTTTG")))   
#unassigned <- as.data.frame(WhichCells(subset, idents = c("NA")))
perf1x$misfit <- "perf1x"
names(perf1x)[1] <- "barcode"
perf2x$misfit <- "perf2x"
names(perf2x)[1] <- "barcode"
#perf4x$misfit <- "perf4x"
#names(perf4x)[1] <- "barcode"
scr$misfit <- "scr"
names(scr)[1] <- "barcode"
empty$misfit <- "empty"
names(empty)[1] <- "barcode"
v1$misfit <- "v1"
names(v1)[1] <- "barcode"
v2$misfit <- "v2"
names(v2)[1] <- "barcode"
v3$misfit <- "v3"
names(v3)[1] <- "barcode"
v4$misfit <- "v4"
names(v4)[1] <- "barcode"
v5$misfit <- "v5"
names(v5)[1] <- "barcode"
v6$misfit <- "v6"
names(v6)[1] <- "barcode"
v7$misfit <- "v7"
names(v7)[1] <- "barcode"
v8$misfit <- "v8"
names(v8)[1] <- "barcode"
v9$misfit <- "v9"
names(v9)[1] <- "barcode"
v10$misfit <- "v10"
names(v10)[1] <- "barcode"
v11$misfit <- "v11"
names(v11)[1] <- "barcode"
v12$misfit <- "v12"
names(v12)[1] <- "barcode"
v13$misfit <- "v13"
names(v13)[1] <- "barcode"
v14$misfit <- "v14"
names(v14)[1] <- "barcode"
v15$misfit <- "v15"
names(v15)[1] <- "barcode"

misfits <- dplyr::bind_rows(perf1x, perf2x)
#misfits <- dplyr::bind_rows(misfits, perf4x) #Due to 0 cells being annotated with 4x Perfect, discard
misfits <- dplyr::bind_rows(misfits, empty)
misfits <- dplyr::bind_rows(misfits, scr)
misfits <- dplyr::bind_rows(misfits, v1)
misfits <- dplyr::bind_rows(misfits, v2)
misfits <- dplyr::bind_rows(misfits, v3)
misfits <- dplyr::bind_rows(misfits, v4)
misfits <- dplyr::bind_rows(misfits, v5)
misfits <- dplyr::bind_rows(misfits, v6)
misfits <- dplyr::bind_rows(misfits, v7)
misfits <- dplyr::bind_rows(misfits, v8)
misfits <- dplyr::bind_rows(misfits, v9)
misfits <- dplyr::bind_rows(misfits, v10)
misfits <- dplyr::bind_rows(misfits, v11)
misfits <- dplyr::bind_rows(misfits, v12)
misfits <- dplyr::bind_rows(misfits, v13)
misfits <- dplyr::bind_rows(misfits, v14)
misfits <- dplyr::bind_rows(misfits, v15)

#ok let's add this as another column to exp1, our df with the misfit binning
misfclusters <- dplyr::inner_join(misfits, clusters, by = "barcode")
misfcls <- data.frame(misfclusters$misfit, misfclusters$cluster)
names(misfcls)[1] <- "miSFITs"
names(misfcls)[2] <- "Cluster"
#mfcs <- misfcls %>% group_by_all %>% dplyr::count
mfcs <- misfcls %>% dplyr::count(miSFITs, Cluster, sort=TRUE)
names(mfcs)[3] <- "Count"

#Add info to metadata - same as miSFITs metadata but more accessible
myBarcode_sub = rownames(subset@meta.data)
misfitsv = misfits[match(myBarcode_sub,misfits$barcode),]
subset$misfits = misfitsv$misfit
subset@active.ident <- factor(subset@active.ident, levels=c("empty", "perf2x", "perf1x", "v10", "v5", 
                                                            "v3", "v8", "v7", "v9", "v6", "v4", "v15", 
                                                            "v11", "v14", "v2", "v12", "v1", "v13", "scr", "NA")) 
subset$misfits <- factor(subset$misfits,
                         levels = c("empty", "perf2x", "perf1x", "v10", "v5", "v3", "v8", "v7", "v9", 
                                    "v6", "v4", "v15", "v11", "v14", "v2", "v12", "v1", "v13", "scr", "NA"))

Idents(subset) <- "Clustering"
subset@active.ident <- factor(subset@active.ident, levels=c("4","1","8","3","7","6","5","0","2")) 
subset$Clustering <- factor(subset$Clustering,
                            levels = c("4", "1", "8", "3", "7", "6", "5", "0", "2"))

subsetPR <- subset(subset, subset = D52NCAR>0)

###Binning by Percentile and Decile Rank
subsetPR <- subset(subset, subset = D52NCAR>0)
carexp <- FetchData(subsetPR, vars =c("D52NCAR"))
carexp <- tibble::rownames_to_column(carexp, "barcode")
car <- carexp
library(dplyr)
car10 = mutate(car, decile_rank = ntile(car$D52NCAR,10))
car100 = mutate(car, percentile_rank = ntile(car$D52NCAR,100))
#assigning decile and percentile rank with clusters 
decile <- dplyr::inner_join(car10, clusters, by = "barcode")
percentile <- dplyr::inner_join(car100, clusters, by = "barcode")

#decile data wrangling
d <- data.frame(decile$decile_rank,decile$cluster)
names(d)[1] <- "Decile"
names(d)[2] <- "Cluster"
decile <- d %>% dplyr::count(Decile, Cluster, sort=TRUE)
#decile <- d %>% group_by_all %>% dplyr::count
names(decile)[3] <- "Count"
#percentile data wrangling
p1 <- dplyr::filter(percentile, percentile_rank >= 1 & percentile_rank <= 10)
p2 <- dplyr::filter(percentile, percentile_rank > 10 & percentile_rank <= 20)
p3 <- dplyr::filter(percentile, percentile_rank > 20 & percentile_rank <= 30)
p4 <- dplyr::filter(percentile, percentile_rank > 30 & percentile_rank <= 40)
p5 <- dplyr::filter(percentile, percentile_rank > 40 & percentile_rank <= 50)
p6 <- dplyr::filter(percentile, percentile_rank > 50 & percentile_rank <= 60)
p7 <- dplyr::filter(percentile, percentile_rank > 60 & percentile_rank <= 70)
p8 <- dplyr::filter(percentile, percentile_rank > 70 & percentile_rank <= 80)
p9 <- dplyr::filter(percentile, percentile_rank > 80 & percentile_rank <= 90)
p10 <- dplyr::filter(percentile, percentile_rank > 90 & percentile_rank <= 100)

p1$p_group <- "1-10"
p2$p_group <- "11-20"
p3$p_group <- "21-30"
p4$p_group <- "31-40"
p5$p_group <- "41-50"
p6$p_group <- "51-60"
p7$p_group <- "61-70"
p8$p_group <- "71-80"
p9$p_group <- "81-90"
p10$p_group <- "91-100"

prank <- dplyr::bind_rows(p1, p2)
prank <- dplyr::bind_rows(prank, p3)
prank <- dplyr::bind_rows(prank, p4)
prank <- dplyr::bind_rows(prank, p5)
prank <- dplyr::bind_rows(prank, p6)
prank <- dplyr::bind_rows(prank, p7)
prank <- dplyr::bind_rows(prank, p8)
prank <- dplyr::bind_rows(prank, p9)
prank <- dplyr::bind_rows(prank, p10)

#we can add this to the metadata
perc_rank <- prank[-c(2:4)]
names(perc_rank)[2] <- "Percentile_Rank"

myBarcode_sub = rownames(subsetPR@meta.data)
Percentile_Rank = perc_rank[match(myBarcode_sub,perc_rank$barcode),]
subsetPR$Percentile_Rank = Percentile_Rank$Percentile_Rank
Idents(subsetPR) <- "Percentile_Rank"
subsetPR@active.ident <- factor(subset@active.ident, levels=c("1-10", "11-20", "21-30", "31-40", "41-50", 
                                                              "51-60", "61-70", "71-80", "81-90", "91-100")) 
subsetPR$Percentile_Rank <- factor(subsetPR$Percentile_Rank,
                                   levels = c("1-10", "11-20", "21-30", "31-40", "41-50", 
                                              "51-60", "61-70", "71-80", "81-90", "91-100"))


###Visualisations 

#Genetic demultiplexing
#As filtered out "DOUBLETS" (e.g., doublets/multiplets and empty droplets) from subset, 
#create subset seuobj with DOUBLETS in it with necessary analysis:   
{
seuobj <- subset(seuobj, 
                 subset = nFeature_RNA >200 & 
                   nFeature_RNA <6000 & percent.mt <10 & 
                   (ids=="A" | ids=="B" | ids=="C" | ids == "DOUBLET"))
seuobj <- NormalizeData(seuobj, normalization.method = "LogNormalize", scale.factor = 10000)
#id variable features
seuobj <- FindVariableFeatures(seuobj, selection.method = "vst")
#view variable features
top50sub <- head(VariableFeatures(seuobj), 50)
#scale data
all.genes.seuobj <- rownames(seuobj)
seuobj <- ScaleData(seuobj, features=all.genes.seuobj)
#PCA
seuobj <- RunPCA(seuobj, features = VariableFeatures(object = seuobj))
seuobj <- FindNeighbors(seuobj, dims=1:50)
seuobj <- FindClusters(seuobj, resolution = 0.5)

#UMAP & TSNE
seuobj <- RunUMAP(seuobj, dims = 1:50)
seuobj <- RunTSNE(seuobj, dims = 1:50)
#Adding to metadata
seuobj[["Clustering"]] <-Idents(object = seuobj)
Idents(seuobj) <- "ids"
}
#Extended Data Figure 11A
Idents(seuobj) <- "ids"
UMAPPlot(seuobj) #Note - DOUBLETS = filtered (e.g., doublets/multiplets & empty droplets)

#Extended Data Figure 11B
Idents(subset) <- "ids"
UMAPPlot(subset)

ANo <- length(WhichCells(seuobj, idents =   c("A")))
BNo <- length(WhichCells(seuobj, idents =   c("B")))
CNo <- length(WhichCells(seuobj, idents =   c("C")))
Filtered <- length(WhichCells(seuobj, idents =   c("DOUBLET")))

GDStatus <- c("Donor A", "Donor B", "Donor C", "Filtered")
GDNumbers <- c(ANo, BNo, CNo, Filtered)

GDRes <- data.frame(GDStatus, GDNumbers)

#Extended Data Figure 11C
ggplot(GDRes, aes(x=GDStatus, y=GDNumbers, fill = GDStatus)) + 
  geom_bar(stat = "identity", color = "black") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none" ) + 
  labs(y= "Number of cells", 
       title = "Genetic Demultiplexing Results") +
  scale_fill_manual(values = c("#00BA38",
                               "#619CFF",
                               "#F8766D",
                               "#C77CFF")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank()) +
  geom_text(aes(label = GDNumbers), size = 3, position = position_stack(vjust = 0.5))

#Extended Data Figure 12F
subsetmisf <- subset(subset, 
                 subset =  misfits=="empty" | misfits=="perf2x" | misfits=="v10" | misfits=="perf1x" | misfits=="v5" | misfits=="v3" |
                              misfits=="v8" | misfits=="v7" | misfits=="v9" | misfits=="v6" | misfits=="v4" | misfits=="v15"| 
                              misfits=="v11" | misfits=="v14" | misfits=="v2" | misfits=="v12" | misfits=="v1" | misfits=="scr"| 
                              misfits=="v13" )

UMAPPlot(subsetmisf, cols = c("empty" =  "#0B0405FF",
                              "perf2x" = "#1A0E19FF",
                              "v10" = "#28192FFF",
                              "perf1x" = "#342346FF",
                              "v5" = "#3B2F5EFF",
                              "v3" = "#403B78FF",
                              "v8" = "#40498EFF",
                              "v7" = "#3A599AFF",
                              "v9" = "#366A9FFF",
                              "v6" = "#357BA2FF",
                              "v4" = "#348AA6FF",
                              "v15" = "#349AAAFF",
                              "v11" = "#38AAACFF",
                              "v14" = "#42B9ADFF",
                              "v2" = "#54C9ADFF",
                              "v12" = "#78D6AEFF",
                              "v1" = "#A0DFB9FF",
                              "scr" = "#79b285",
                              "v13" = "#DEF5E5FF"))

##Extended Data Figure 12G
#Left
VlnPlot(subsetmisf, features = "D52NCAR", pt.size = 0 , cols = c("empty" =  "#0B0405FF",
                                                   "perf2x" = "#1A0E19FF",
                                                   "v10" = "#28192FFF",
                                                   "perf1x" = "#342346FF",
                                                   "v5" = "#3B2F5EFF",
                                                   "v3" = "#403B78FF",
                                                   "v8" = "#40498EFF",
                                                   "v7" = "#3A599AFF",
                                                   "v9" = "#366A9FFF",
                                                   "v6" = "#357BA2FF",
                                                   "v4" = "#348AA6FF",
                                                   "v15" = "#349AAAFF",
                                                   "v11" = "#38AAACFF",
                                                   "v14" = "#42B9ADFF",
                                                   "v2" = "#54C9ADFF",
                                                   "v12" = "#78D6AEFF",
                                                   "v1" = "#A0DFB9FF",
                                                   "scr" = "#79b285",
                                                   "v13" = "#DEF5E5FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0) +  
  NoLegend()

#Right 
CARexpression <- FetchData(subset, vars =c("D52NCAR"))
CARexpression <- tibble::rownames_to_column(CARexpression, "barcode")

CAR_perf1x <- dplyr::semi_join(CARexpression, perf1x, by = "barcode")
CAR_perf2x <- dplyr::semi_join(CARexpression, perf2x, by = "barcode")
#CAR_perf4x <- dplyr::semi_join(CARexpression, perf4x, by = "barcode")
CAR_scr <- dplyr::semi_join(CARexpression, scr, by = "barcode")
CAR_empty <- dplyr::semi_join(CARexpression, empty, by = "barcode")
CAR_v1 <- dplyr::semi_join(CARexpression, v1, by = "barcode")
CAR_v2 <- dplyr::semi_join(CARexpression, v2, by = "barcode")
CAR_v3 <- dplyr::semi_join(CARexpression, v3, by = "barcode")
CAR_v4 <- dplyr::semi_join(CARexpression, v4, by = "barcode")
CAR_v5 <- dplyr::semi_join(CARexpression, v5, by = "barcode")
CAR_v6 <- dplyr::semi_join(CARexpression, v6, by = "barcode")
CAR_v7 <- dplyr::semi_join(CARexpression, v7, by = "barcode")
CAR_v8 <- dplyr::semi_join(CARexpression, v8, by = "barcode")
CAR_v9 <- dplyr::semi_join(CARexpression, v9, by = "barcode")
CAR_v10 <- dplyr::semi_join(CARexpression, v10, by = "barcode")
CAR_v11 <- dplyr::semi_join(CARexpression, v11, by = "barcode")
CAR_v12 <- dplyr::semi_join(CARexpression, v12, by = "barcode")
CAR_v13 <- dplyr::semi_join(CARexpression, v13, by = "barcode")
CAR_v14 <- dplyr::semi_join(CARexpression, v14, by = "barcode")
CAR_v15 <- dplyr::semi_join(CARexpression, v15, by = "barcode")

CAR_perf1x$miSFIT <- "Perfect1x"
CAR_perf2x$miSFIT <- "Perfect2x"
CAR_empty$miSFIT <- "Empty"
CAR_scr$miSFIT <- "Scramble"
CAR_v1$miSFIT <- "V1"
CAR_v2$miSFIT <- "V2"
CAR_v3$miSFIT <- "V3"
CAR_v4$miSFIT <- "V4"
CAR_v5$miSFIT <- "V5"
CAR_v6$miSFIT <- "V6"
CAR_v7$miSFIT <- "V7"
CAR_v8$miSFIT <- "V8"
CAR_v9$miSFIT <- "V9"
CAR_v10$miSFIT <- "V10"
CAR_v11$miSFIT <- "V11"
CAR_v12$miSFIT <- "V12"
CAR_v13$miSFIT <- "V13"
CAR_v14$miSFIT <- "V14"
CAR_v15$miSFIT <- "V15"

CARperf1x_mean <- CAR_perf1x %>% summarise_if(is.numeric, mean)
CARperf2x_mean <- CAR_perf2x %>% summarise_if(is.numeric, mean)
#CARperf4x_mean <- CAR_perf4x %>% summarise_if(is.numeric, mean)
CARscr_mean <- CAR_scr %>% summarise_if(is.numeric, mean)
CARempty_mean <- CAR_empty %>% summarise_if(is.numeric, mean)
CARv1_mean <- CAR_v1 %>% summarise_if(is.numeric, mean)
CARv2_mean <- CAR_v2 %>% summarise_if(is.numeric, mean)
CARv3_mean <- CAR_v3 %>% summarise_if(is.numeric, mean)
CARv4_mean <- CAR_v4 %>% summarise_if(is.numeric, mean)
CARv5_mean <- CAR_v5 %>% summarise_if(is.numeric, mean)
CARv6_mean <- CAR_v6 %>% summarise_if(is.numeric, mean)
CARv7_mean <- CAR_v7 %>% summarise_if(is.numeric, mean)
CARv8_mean <- CAR_v8 %>% summarise_if(is.numeric, mean)
CARv9_mean <- CAR_v9 %>% summarise_if(is.numeric, mean)
CARv10_mean <- CAR_v10 %>% summarise_if(is.numeric, mean)
CARv11_mean <- CAR_v11 %>% summarise_if(is.numeric, mean)
CARv12_mean <- CAR_v12 %>% summarise_if(is.numeric, mean)
CARv13_mean <- CAR_v13 %>% summarise_if(is.numeric, mean)
CARv14_mean <- CAR_v14 %>% summarise_if(is.numeric, mean)
CARv15_mean <- CAR_v15 %>% summarise_if(is.numeric, mean)

rownames(CARperf1x_mean) <- c("Perfect1x")
rownames(CARperf2x_mean) <- c("Perfect2x")
#rownames(CARperf4x_mean) <- c("Perfect4x")
rownames(CARscr_mean) <- c("Scramble")
rownames(CARempty_mean) <- c("Empty")
rownames(CARv1_mean) <- c("V1")
rownames(CARv2_mean) <- c("V2")
rownames(CARv3_mean) <- c("V3")
rownames(CARv4_mean) <- c("V4")
rownames(CARv5_mean) <- c("V5")
rownames(CARv6_mean) <- c("V6")
rownames(CARv7_mean) <- c("V7")
rownames(CARv8_mean) <- c("V8")
rownames(CARv9_mean) <- c("V9")
rownames(CARv10_mean) <- c("V10")
rownames(CARv11_mean) <- c("V11")
rownames(CARv12_mean) <- c("V12")
rownames(CARv13_mean) <- c("V13")
rownames(CARv14_mean) <- c("V14")
rownames(CARv15_mean) <- c("V15")

CAR_means <- dplyr::bind_rows(CARperf1x_mean, CARperf2x_mean)
#CAR_means <- dplyr::bind_rows(CAR_means, CARperf4x_mean) 
CAR_means <- dplyr::bind_rows(CAR_means, CARscr_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARempty_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv1_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv2_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv3_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv4_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv5_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv6_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv7_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv8_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv9_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv10_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv11_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv12_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv13_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv14_mean)
CAR_means <- dplyr::bind_rows(CAR_means, CARv15_mean)
CAR_means <- tibble::rownames_to_column(CAR_means, "miSFIT")

CAR_means <- CAR_means[order( CAR_means[,2] ),]
CARm_order <- CAR_means$miSFIT
CAR_means$D52NCAR <- signif(CAR_means$D52NCAR,3)

ggplot(CAR_means, aes(x = factor(miSFIT, level = CARm_order), y= D52NCAR, fill = miSFIT)) + 
geom_bar(stat = "identity",colour = "black") + 
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank()) + 
labs(x= "miSFIT", y= "Average CAR Expression", title= "miSFIT CAR Expression") +
geom_text(aes(label=D52NCAR), vjust=-0.15, color="black", size=3.5) +
scale_fill_manual( values = c("Empty" =  "#0B0405FF",
                              "Perfect2x" = "#1A0E19FF",
                              "V10" = "#28192FFF",
                              "Perfect1x" = "#342346FF",
                              "V5" = "#3B2F5EFF",
                              "V3" = "#403B78FF",
                              "V8" = "#40498EFF",
                              "V7" = "#3A599AFF",
                              "V9" = "#366A9FFF",
                              "V6" = "#357BA2FF",
                              "V4" = "#348AA6FF",
                              "V15" = "#349AAAFF",
                              "V11" = "#38AAACFF",
                              "V14" = "#42B9ADFF",
                              "V2" = "#54C9ADFF",
                              "V12" = "#78D6AEFF",
                              "V1" = "#A0DFB9FF",
                              "Scramble" = "#79b285",
                              "V13" = "#DEF5E5FF")) + NoLegend()

#Figure 3B
FeaturePlot(subset, features = "D52NCAR")

Idents(subsetPR) <- "Percentile_Rank"
#Figure 3C
UMAPPlot(subsetPR, cols = c( "#440154FF",
                             "#482878FF",
                             "#3E4A89FF",
                             "#31688EFF",
                             "#26828EFF",
                             "#1F9E89FF",
                             "#35B779FF",
                             "#6DCD59FF",
                             "#B4DE2CFF",
                             "#FDE725FF"))


##Figure 3D
#Left
RidgePlot(subsetPR, log = TRUE, features = "D52NCAR",cols = c( "#440154FF",
                                                   "#482878FF",
                                                   "#3E4A89FF",
                                                   "#31688EFF",
                                                   "#26828EFF",
                                                   "#1F9E89FF",
                                                   "#35B779FF",
                                                   "#6DCD59FF",
                                                   "#B4DE2CFF",
                                                   "#FDE725FF"))
#Right
CARp1 <- dplyr::semi_join(carexp,p1, by = "barcode")
CARp2 <- dplyr::semi_join(carexp,p2, by = "barcode")
CARp3 <- dplyr::semi_join(carexp,p3, by = "barcode")
CARp4 <- dplyr::semi_join(carexp,p4, by = "barcode")
CARp5 <- dplyr::semi_join(carexp,p5, by = "barcode")
CARp6 <- dplyr::semi_join(carexp,p6, by = "barcode")
CARp7 <- dplyr::semi_join(carexp,p7, by = "barcode")
CARp8 <- dplyr::semi_join(carexp,p8, by = "barcode")
CARp9 <- dplyr::semi_join(carexp,p9, by = "barcode")
CARp10 <- dplyr::semi_join(carexp,p10, by = "barcode")

CARp1m <- CARp1 %>% summarise_if(is.numeric, mean)
CARp2m <- CARp2 %>% summarise_if(is.numeric, mean)
CARp3m <- CARp3 %>% summarise_if(is.numeric, mean)
CARp4m <- CARp4 %>% summarise_if(is.numeric, mean)
CARp5m <- CARp5 %>% summarise_if(is.numeric, mean)
CARp6m <- CARp6 %>% summarise_if(is.numeric, mean)
CARp7m <- CARp7 %>% summarise_if(is.numeric, mean)
CARp8m <- CARp8 %>% summarise_if(is.numeric, mean)
CARp9m <- CARp9 %>% summarise_if(is.numeric, mean)
CARp10m <- CARp10 %>% summarise_if(is.numeric, mean)

CARpercentile <- rbind(CARp1m,CARp2m, CARp3m,CARp4m,CARp5m,CARp6m,
                       CARp7m,CARp8m,CARp9m, CARp10m)
CARpercentile$Percentile <- c("Percentile 1-10", "Percentile 11-20", "Percentile 21-30","Percentile 31-40",
                              "Percentile 41-50", "Percentile 51-60","Percentile 61-70","Percentile 71-80",
                              "Percentile 81-90","Percentile 91-100")

CARp1$Percentile <- "Percentile 1-10"
CARp2$Percentile <- "Percentile 11-20"
CARp3$Percentile <- "Percentile 21-30"
CARp4$Percentile <- "Percentile 31-40"
CARp5$Percentile <- "Percentile 41-50"
CARp6$Percentile <- "Percentile 51-60"
CARp7$Percentile <- "Percentile 61-70"
CARp8$Percentile <- "Percentile 71-80"
CARp9$Percentile <- "Percentile 81-90"
CARp10$Percentile <- "Percentile 91-100"

percentile_total <- sum(length(p1$barcode), length(p2$barcode), length(p3$barcode), length(p4$barcode),
                        length(p5$barcode), length(p6$barcode), length(p7$barcode), length(p8$barcode),
                        length(p9$barcode), length(p10$barcode))
percentilebin_average <- mean(length(p1$barcode), length(p2$barcode), length(p3$barcode), length(p4$barcode),
                              length(p5$barcode), length(p6$barcode), length(p7$barcode), length(p8$barcode),
                              length(p9$barcode), length(p10$barcode))

ggplot(CARpercentile, aes(x = Percentile, y = D52NCAR, fill = Percentile)) + 
  geom_bar(stat = "identity") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none" ) + 
  labs(x= "Percentile Rank Bin", y= "Average CAR Expression") +
  scale_fill_manual(values = c("Percentile 1-10" =  "#440154FF",
                               "Percentile 11-20" = "#482878FF",
                               "Percentile 21-30" = "#3E4A89FF",
                               "Percentile 31-40" = "#31688EFF",
                               "Percentile 41-50" = "#26828EFF",
                               "Percentile 51-60" = "#1F9E89FF",
                               "Percentile 61-70" = "#35B779FF",
                               "Percentile 71-80" = "#6DCD59FF",
                               "Percentile 81-90" = "#B4DE2CFF",
                               "Percentile 91-100" = "#FDE725FF"))

#Figure 3e
VlnPlot(subsetPR, features = c("CD69"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                    "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                    "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                    "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0) +  
  NoLegend() +
  labs(title = "CD69")

VlnPlot(subsetPR, features = c("CCL5"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0) +  
  NoLegend() +
  labs(title = "CCL5")

VlnPlot(subsetPR, features = c("LAG3"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0) +  
  NoLegend() +
  labs(title = "LAG3")

VlnPlot(subsetPR, features = c("CCR7"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0) +  
  NoLegend() +
  labs(title = "CCR7")

