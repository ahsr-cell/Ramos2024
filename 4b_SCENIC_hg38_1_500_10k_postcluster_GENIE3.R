###SCENIC: Post GENIE3 steps - to completion of SCENIC analysis

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

setwd("/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/hg38/CCB_25JulyRun/refseq")
knitr::opts_knit$set(root.dir="/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/hg38/CCB_25JulyRun/refseq") 

scenicOptions <- readRDS("/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/hg38/CCB_25JulyRun/refseq/int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 8
scenicOptions@settings$seed <- 123

load("/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/18July_hg38/PreProcessing_for_hg38SCENIC_18July.RData")
exprMat <- as.matrix(exprMat)
#exprMat <- readRDS("/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/18July_hg38/exprMat.rds")
#cellInfo <- data.frame(subsetPR$Percentile_Rank, stringsAsFactors = FALSE)
#scenicOptions <- initializeScenic(org="hgnc",
#                                  dbDir = "/Users/ahsramos/Desktop/scRNAseq_bioinformatic_analysis/SCENIC/cisTarget_databases/hg38",
#                                  dbs=c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", 
#                                        "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"),
#                                  datasetTitle="D52NCAR_PR_hg38_100-500bp_and_10kb", nCores=8) 
#scenicOptions@settings$db_mcVersion <- "v9"
#scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
#saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #done
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) #done
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat) #done
 
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat) #done
savedSelections <- shiny::runApp(aucellApp) #done
newThresholds <- savedSelections$thresholds #done
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds" #done
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds")) #done
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat = exprMat) #done

Regulon_Rankings <- loadInt(scenicOptions,"aucell_rankings")
RegulonAUC <- readRDS("int/3.4_regulonAUC.Rds")
RegulonAUCThresholds <- readRDS("int/3.5_AUCellThresholds.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byPR <- sapply(split(rownames(cellInfo), cellInfo$Percentile_Rank),
                               function(cells) rowMeans(getAUC(regulonAUC)[,cells]))


###Filter identified regulons using Mutual Information against CAR expression 
##Necessary data wrangling
regAUC <- regulonAUC@assays@data@listData[["AUC"]] 

regulonAUC_values <- as.data.frame(t(regAUC), stringsAsFactors = FALSE)
regulonAUC_values <- tibble::rownames_to_column(regulonAUC_values, "barcode")

UBRegAUCPR1 <- dplyr::semi_join(regulonAUC_values, p1, by = "barcode")
UBRegAUCPR2 <- dplyr::semi_join(regulonAUC_values, p2, by = "barcode")
UBRegAUCPR3 <- dplyr::semi_join(regulonAUC_values, p3, by = "barcode")
UBRegAUCPR4 <- dplyr::semi_join(regulonAUC_values, p4, by = "barcode")
UBRegAUCPR5 <- dplyr::semi_join(regulonAUC_values, p5, by = "barcode")
UBRegAUCPR6 <- dplyr::semi_join(regulonAUC_values, p6, by = "barcode")
UBRegAUCPR7 <- dplyr::semi_join(regulonAUC_values, p7, by = "barcode")
UBRegAUCPR8 <- dplyr::semi_join(regulonAUC_values, p8, by = "barcode")
UBRegAUCPR9 <- dplyr::semi_join(regulonAUC_values, p9, by = "barcode")
UBRegAUCPR10 <- dplyr::semi_join(regulonAUC_values, p10, by = "barcode")

UBRegAUCPR1_mean <- UBRegAUCPR1 %>% summarise_if(is.numeric, mean)
UBRegAUCPR2_mean <- UBRegAUCPR2 %>% summarise_if(is.numeric, mean)
UBRegAUCPR3_mean <- UBRegAUCPR3 %>% summarise_if(is.numeric, mean)
UBRegAUCPR4_mean <- UBRegAUCPR4 %>% summarise_if(is.numeric, mean)
UBRegAUCPR5_mean <- UBRegAUCPR5 %>% summarise_if(is.numeric, mean)
UBRegAUCPR6_mean <- UBRegAUCPR6 %>% summarise_if(is.numeric, mean)
UBRegAUCPR7_mean <- UBRegAUCPR7 %>% summarise_if(is.numeric, mean)
UBRegAUCPR8_mean <- UBRegAUCPR8 %>% summarise_if(is.numeric, mean)
UBRegAUCPR9_mean <- UBRegAUCPR9 %>% summarise_if(is.numeric, mean)
UBRegAUCPR10_mean <- UBRegAUCPR10 %>% summarise_if(is.numeric, mean)

UBRegAUCPR1_mean <- as.data.frame(UBRegAUCPR1_mean, stringsAsFactors = FALSE)
UBRegAUCPR2_mean <- as.data.frame(UBRegAUCPR2_mean, stringsAsFactors = FALSE)
UBRegAUCPR3_mean <- as.data.frame(UBRegAUCPR3_mean, stringsAsFactors = FALSE)
UBRegAUCPR4_mean <- as.data.frame(UBRegAUCPR4_mean, stringsAsFactors = FALSE)
UBRegAUCPR5_mean <- as.data.frame(UBRegAUCPR5_mean, stringsAsFactors = FALSE)
UBRegAUCPR6_mean <- as.data.frame(UBRegAUCPR6_mean, stringsAsFactors = FALSE)
UBRegAUCPR7_mean <- as.data.frame(UBRegAUCPR7_mean, stringsAsFactors = FALSE)
UBRegAUCPR8_mean <- as.data.frame(UBRegAUCPR8_mean, stringsAsFactors = FALSE)
UBRegAUCPR9_mean <- as.data.frame(UBRegAUCPR9_mean, stringsAsFactors = FALSE)
UBRegAUCPR10_mean <- as.data.frame(UBRegAUCPR10_mean, stringsAsFactors = FALSE)

UBPRMRegAUCPR_MI <- rbind.data.frame(UBRegAUCPR1_mean, UBRegAUCPR2_mean, UBRegAUCPR3_mean, UBRegAUCPR4_mean, UBRegAUCPR5_mean,
                                     UBRegAUCPR6_mean, UBRegAUCPR7_mean, UBRegAUCPR8_mean, UBRegAUCPR9_mean, UBRegAUCPR10_mean, stringsAsFactors = FALSE)

UBPRMRegAUCPR_MI <- cbind.data.frame(UBPRMRegAUCPR_MI, t(UBPRMGenesa[which(rownames(UBPRMGenesa) == "D52NCAR"),]), stringsAsFactors = FALSE)

##Applying calculating Mutual Information
require(infotheo)    
discret_UBPRMRegAUCPR_MI <- discretize(UBPRMRegAUCPR_MI)
rownames(discret_UBPRMRegAUCPR_MI) <- rownames(UBPRMRegAUCPR_MI)
CARdisc_Reg <- discretize(UBPRMRegAUCPR_MI$D52NCAR)
RegCARE <- infotheo::entropy(CARdisc_Reg)
calc_NMI <- NULL
for (i in 1:(which(colnames(UBPRMRegAUCPR_MI) == "D52NCAR"))) {
  OIdisc_Reg <- discretize(UBPRMRegAUCPR_MI[,i])
  RegOIE <- infotheo::entropy(OIdisc_Reg)
  MI <- mutinformation(CARdisc_Reg, OIdisc_Reg)
  NMI <- (MI / sqrt(RegCARE*RegOIE))
  a <- cbind.data.frame(MI,NMI)
  calc_NMI <- rbind.data.frame(calc_NMI, a)
}
calc_RegNMI <- calc_NMI
rownames(calc_RegNMI) <- colnames(UBPRMRegAUCPR_MI)
calc_RegNMI1 <- dplyr::filter(calc_RegNMI, NMI == "1") #Filtering for highest mutual information
calc_RegNMI1 <- rownames_to_column(calc_RegNMI1, "Regulon")

#Figure 4e 
regulonActivity_byPR_Scaled <- t(scale(t(regulonActivity_byPR), center = T, scale=T))
regulonActivity_byPR_Scaled_1NMI <- as.data.frame(regulonActivity_byPR_Scaled, stringsAsFactors = F)
regulonActivity_byPR_Scaled_1NMI <- rownames_to_column(data.frame(regulonActivity_byPR_Scaled), "Regulon")
regulonActivity_byPR_Scaled_1NMI <- semi_join(regulonActivity_byPR_Scaled_1NMI, calc_RegNMI1, "Regulon")
colnames(regulonActivity_byPR_Scaled_1NMI)[1] <- "Regulon"
regulonActivity_byPR_Scaled_1NMI <- column_to_rownames(regulonActivity_byPR_Scaled_1NMI, "Regulon")
colnames(regulonActivity_byPR_Scaled_1NMI) <- c("1-10", "11-20", "21-30", "31-40", "41-50",
                                            "51-60", "61-70", "71-80", "81-90", "91-100")
ComplexHeatmap::Heatmap(as.matrix(regulonActivity_byPR_Scaled_1NMI), name="Regulon activity",
                        column_order = c("1-10", "11-20", "21-30", "31-40", "41-50",
                                         "51-60", "61-70", "71-80", "81-90", "91-100"))

###Extended Data Figure 17
regulonActivity_byPR <- rownames_to_column(data.frame(regulonActivity_byPR), "Regulon")
regulonActivity_byPR_reg1NMI_linear <-  semi_join(regulonActivity_byPR, calc_RegNMI1, "Regulon")
regulonActivity_byPR_reg1NMI_linear <- column_to_rownames(regulonActivity_byPR_reg1NMI_linear, "Regulon")
colnames(regulonActivity_byPR_reg1NMI_linear) <- c("1-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
holder <- UBPRMGenesa[which(rownames(UBPRMGenesa) == "D52NCAR"),]
colnames(holder) <- c("1-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")

regulonActivity_byPR_reg1NMI_linear <- rbind.data.frame(regulonActivity_byPR_reg1NMI_linear, holder, stringsAsFactors = FALSE)

car_index <- which(rownames(regulonActivity_byPR_reg1NMI_linear) == "D52NCAR")  #CAR is the last row
length_sc <- length(rownames(regulonActivity_byPR_reg1NMI_linear))-1 #subtracted one as our CAR is the last row
unbiased_df_lin_PR_reg <- NULL                 #initialise empty df to fill                                       
for (i in 1:length_sc) {
  ub_linfit_PR <- lm(unlist(regulonActivity_byPR_reg1NMI_linear[i,]) ~ unlist(regulonActivity_byPR_reg1NMI_linear[car_index,])) #conduct linear fit
  goi_PR <- rownames(regulonActivity_byPR_reg1NMI_linear[i,]) #gene name
  r2_PR <- as.numeric(summary(ub_linfit_PR)$r.squared) #Rsquared
  pval_PR <- as.numeric(summary(ub_linfit_PR)$coefficients[2,4]) #p value
  linf_PR_reg <- cbind(goi_PR, r2_PR, pval_PR) #create new row
  unbiased_df_lin_PR_reg <- rbind(unbiased_df_lin_PR_reg, linf_PR_reg) #add to initialised df 
}                                 
unbiased_df_lin_PR_reg <- as.data.frame(unbiased_df_lin_PR_reg, stringsAsFactors=FALSE) #create intermediate df
udlPR_goi <- unbiased_df_lin_PR_reg$goi
udlPR_r2 <- unbiased_df_lin_PR_reg$r2 #convert from character to numeric
udlPR_r2 <- as.numeric(udlPR_r2)
udlPR_p <- unbiased_df_lin_PR_reg$pval
udlPR_p <- as.numeric(udlPR_p) #convert from character to numeric
unbiased_df_lin_PR_reg <- data.frame(udlPR_goi, udlPR_r2, udlPR_p, stringsAsFactors = FALSE) #create final df
names(unbiased_df_lin_PR_reg)[1] <- "Regulon"
names(unbiased_df_lin_PR_reg)[2] <- "Rsquared"
names(unbiased_df_lin_PR_reg)[3] <- "P"
unbiased_df_lin_PR_reg_sansNA <- na.omit(unbiased_df_lin_PR_reg) #drop anything with R2 or P with NaN
length(which(unbiased_df_lin_PR_reg_sansNA$Rsquared>=0.8))

ggplot(unbiased_df_lin_PR_reg_sansNA, aes(x=Rsquared)) +
  theme_minimal() +
  geom_histogram(binwidth=0.03, fill = "steelblue") +
  xlab("R2") +
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

