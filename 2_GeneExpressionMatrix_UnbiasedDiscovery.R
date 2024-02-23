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
require(readxl)
require(escape)
require(dittoSeq)
require(GSVA)
require(biomaRt)
require(EnrichmentBrowser)
require(scales)
require(magrittr)
require(infotheo)
require(wesanderson)

###Conduct 1_Startup[...].R script for necessary seurat object subsetting, data wrangling, etc. 

#Create gene expression matrix based on CAR Expression bins 
UBGenesPR <- FetchData(subsetPR, vars = rownames(subsetPR))
UBGenesPR <- tibble::rownames_to_column(UBGenesPR, "barcode")

UBGenesPR1 <- dplyr::semi_join(UBGenesPR, p1, by = "barcode")
UBGenesPR2 <- dplyr::semi_join(UBGenesPR, p2, by = "barcode")
UBGenesPR3 <- dplyr::semi_join(UBGenesPR, p3, by = "barcode")
UBGenesPR4 <- dplyr::semi_join(UBGenesPR, p4, by = "barcode")
UBGenesPR5 <- dplyr::semi_join(UBGenesPR, p5, by = "barcode")
UBGenesPR6 <- dplyr::semi_join(UBGenesPR, p6, by = "barcode")
UBGenesPR7 <- dplyr::semi_join(UBGenesPR, p7, by = "barcode")
UBGenesPR8 <- dplyr::semi_join(UBGenesPR, p8, by = "barcode")
UBGenesPR9 <- dplyr::semi_join(UBGenesPR, p9, by = "barcode")
UBGenesPR10 <- dplyr::semi_join(UBGenesPR, p10, by = "barcode")

UBGenesPR1_mean <- UBGenesPR1 %>% summarise_if(is.numeric, mean)
UBGenesPR2_mean <- UBGenesPR2 %>% summarise_if(is.numeric, mean)
UBGenesPR3_mean <- UBGenesPR3 %>% summarise_if(is.numeric, mean)
UBGenesPR4_mean <- UBGenesPR4 %>% summarise_if(is.numeric, mean)
UBGenesPR5_mean <- UBGenesPR5 %>% summarise_if(is.numeric, mean)
UBGenesPR6_mean <- UBGenesPR6 %>% summarise_if(is.numeric, mean)
UBGenesPR7_mean <- UBGenesPR7 %>% summarise_if(is.numeric, mean)
UBGenesPR8_mean <- UBGenesPR8 %>% summarise_if(is.numeric, mean)
UBGenesPR9_mean <- UBGenesPR9 %>% summarise_if(is.numeric, mean)
UBGenesPR10_mean <- UBGenesPR10 %>% summarise_if(is.numeric, mean)

rownames(UBGenesPR1_mean) <- c("Percentile_1-10")
rownames(UBGenesPR2_mean) <- c("Percentile_11-20")
rownames(UBGenesPR3_mean) <- c("Percentile_21-30")
rownames(UBGenesPR4_mean) <- c("Percentile_31-40")
rownames(UBGenesPR5_mean) <- c("Percentile_41-50")
rownames(UBGenesPR6_mean) <- c("Percentile_51-60")
rownames(UBGenesPR7_mean) <- c("Percentile_61-70")
rownames(UBGenesPR8_mean) <- c("Percentile_71-80")
rownames(UBGenesPR9_mean) <- c("Percentile_81-90")
rownames(UBGenesPR10_mean) <- c("Percentile_91-100")

UBGenesPR1_meana <- as.data.frame(t(UBGenesPR1_mean), stringsAsFactors = FALSE)
UBGenesPR2_meana <- as.data.frame(t(UBGenesPR2_mean), stringsAsFactors = FALSE)
UBGenesPR3_meana <- as.data.frame(t(UBGenesPR3_mean), stringsAsFactors = FALSE)
UBGenesPR4_meana <- as.data.frame(t(UBGenesPR4_mean), stringsAsFactors = FALSE)
UBGenesPR5_meana <- as.data.frame(t(UBGenesPR5_mean), stringsAsFactors = FALSE)
UBGenesPR6_meana <- as.data.frame(t(UBGenesPR6_mean), stringsAsFactors = FALSE)
UBGenesPR7_meana <- as.data.frame(t(UBGenesPR7_mean), stringsAsFactors = FALSE)
UBGenesPR8_meana <- as.data.frame(t(UBGenesPR8_mean), stringsAsFactors = FALSE)
UBGenesPR9_meana <- as.data.frame(t(UBGenesPR9_mean), stringsAsFactors = FALSE)
UBGenesPR10_meana <- as.data.frame(t(UBGenesPR10_mean), stringsAsFactors = FALSE)

UBPRMGenesa <- cbind(UBGenesPR1_meana, UBGenesPR2_meana, UBGenesPR3_meana, UBGenesPR4_meana, UBGenesPR5_meana,
                     UBGenesPR6_meana, UBGenesPR7_meana, UBGenesPR8_meana, UBGenesPR9_meana, UBGenesPR10_meana)

#Create normalised gene expression matrix by Z-score 
UBPRMGenesa_z <- column_to_rownames(UBPRMGenesa, "GOI")
z_score <- function(x){(x- mean(x, na.rm = FALSE)) /sd(x, na.rm = FALSE)}
UBPW_PRz <- function(UBaveraged) {
  resc_df <- NULL
  rowNumbs <- length(rownames(UBaveraged))
  for (i in 1:rowNumbs){
    require(scales)
    z <- UBPRMGenesa_z %>% dplyr::slice(i)
    resc_row <- z_score(unlist(z))
    resc_df <- rbind(resc_df, resc_row)
  }
  rownames(resc_df) <- rownames(UBaveraged)
  colnames(resc_df) <- colnames(UBaveraged)
  resc_df <- as.data.frame(resc_df, stringsAsFactors = FALSE)
  print("success")
  return(resc_df)
}
UBPRMGenesa_zed <- UBPW_PRz(UBPRMGenesa_z)
UBPRMGenesa_zed <- rownames_to_column(UBPRMGenesa_zed, "GOI")


###Extended Data Figure 13 
Genes <- c("D52NCAR", "CD69", "LAG3","CCR7", "CCL5")
UBPRMGenesa <- column_to_rownames(UBPRMGenesa, "GOI")
UBMod <- UBPRMGenesa[rownames(UBPRMGenesa) %in% Genes, ] 
UBMod <- data.frame(t(UBMod), stringsAsFactors = F)
rownames(UBMod) <- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80", "81-90","91-100")

UBMod1 <- UBPRMGenesa_zed[rownames(UBPRMGenesa_zed) %in% Genes, ] 
UBMod1 <- data.frame(t(UBMod1), stringsAsFactors = F)
rownames(UBMod1) <- c("P1", "P2","P3","P4","P5","P6","P7","P8", "P9","P10")

UBlofa1 <- lm(CD69 ~ log(D52NCAR), data = UBMod)
UBlofa2 <- lm(CCL5 ~ log(D52NCAR), data = UBMod)
UBlfa1 <- lm(LAG3 ~ D52NCAR, data = UBMod)
UBlofa3 <- lm(CCR7 ~ log(D52NCAR), data = UBMod)

#CD69
ggplot(UBMod, aes(x=D52NCAR, y=CD69)) + 
  stat_smooth(method='lm', formula = y ~ log(x),se=F) +
  geom_point(alpha=0.5) + theme_linedraw() +
  geom_text_repel(label=rownames(UBMod), max.overlaps = Inf) +     
  labs(x="", y="",
       title="CD69",
       subtitle= paste("Radj2 =",signif(summary(UBlofa1)$r.squared,4), 
                       "P =", signif(summary(UBlofa1)$coefficients[2,4],4)))
#CCL5
ggplot(UBMod, aes(x=D52NCAR, y=CCL5)) + 
  stat_smooth(method='lm', formula = y ~ log(x),se=F) +
  geom_point(alpha=0.5) + theme_linedraw() +
  geom_text_repel(label=rownames(UBMod), max.overlaps = Inf) +     
  labs(x="", y="",
       title="CCL5",
       subtitle= paste("Radj2 =",signif(summary(UBlofa2)$r.squared,4), 
                       "P =", signif(summary(UBlofa2)$coefficients[2,4],4)))
#LAG3
ggplot(UBMod, aes(x=D52NCAR, y=LAG3)) + 
  stat_smooth(method='lm', formula = y ~ x,se=F) +
  geom_point(alpha=0.5) + theme_linedraw() +
  geom_text_repel(label=rownames(UBMod), max.overlaps = Inf) +
  labs(x="", y="",
       title="LAG3",
       subtitle= paste("Radj2=",signif(summary(UBlfa1)$adj.r.squared,3),
                       ", P =", signif(summary(UBlfa1)$coefficients[2,4],3)))
#CCR7
ggplot(UBMod, aes(x=D52NCAR, y=CCR7)) + 
  stat_smooth(method='lm', formula = y ~ log(x),se=F) +
  geom_point(alpha=0.5) + theme_linedraw() +
  geom_text_repel(label=rownames(UBMod), max.overlaps = Inf) +     
  labs(x="", y="",
          title="CCR7",
       subtitle= paste("Radj2 =",signif(summary(UBlofa3)$r.squared,4), 
                       "P =", signif(summary(UBlofa3)$coefficients[2,4],4)))

###Extended Data Figure 14
#EDF14A
#Import TargetScan for miR17~92
TS_1792 <- read_excel("/Users/ahsramos/Desktop/Manuscript/TargetScan7.2__miR-17-5p_20-5p_93-5p_106-5p_519-3p.predicted_targets.xlsx",
                      sheet = 1)
mi1792t <- TS_1792$`Target gene`
mi1792t <- c(mi1792t, "D52NCAR")
UBPRMGenesa_zed <- rownames_to_column(UBPRMGenesa_zed, "GOI")
mR1792t <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% mi1792t,], stringsAsFactors = F)
rownames(mR1792t) <- mR1792t$GOI
mR1792t <- mR1792t[,-1]
colnames(mR1792t) <- c("1-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
mR1792t <- na.omit(mR1792t)
pheatmap(as.matrix(mR1792t), cluster_cols = F)

#EDF14B
#Conduct linear regression analysis
#linear fitting for miR17~92
{
  car_index <- which(rownames(mR1792t) == "D52NCAR")  #CAR is the last row
  length_sc <- length(rownames(mR1792t))-1 #subtracted one as our CAR is the last row
  unbiased_df_lin_PR <- NULL                 #initialise empty df to fill                                       
  for (i in 1:length_sc) {
    ub_linfit_PR <- lm(unlist(mR1792t[i,]) ~ unlist(mR1792t[car_index,])) #conduct linear fit
    goi_PR <- rownames(mR1792t[i,]) #gene name
    r2_PR <- as.numeric(summary(ub_linfit_PR)$r.squared) #Rsquared
    pval_PR <- as.numeric(summary(ub_linfit_PR)$coefficients[2,4]) #p value
    linf_PR <- cbind(goi_PR, r2_PR, pval_PR) #create new row
    unbiased_df_lin_PR <- rbind(unbiased_df_lin_PR, linf_PR) #add to initialised df 
  }                                 
  unbiased_df_lin_PR <- as.data.frame(unbiased_df_lin_PR, stringsAsFactors=FALSE) #create intermediate df
  udlPR_goi <- unbiased_df_lin_PR$goi
  udlPR_r2 <- unbiased_df_lin_PR$r2 #convert from character to numeric
  udlPR_r2 <- as.numeric(udlPR_r2)
  udlPR_p <- unbiased_df_lin_PR$pval
  udlPR_p <- as.numeric(udlPR_p) #convert from character to numeric
  #summary(M.lm)$adj.r.squared - in case we want to switch to adjustedRsquared
  unbiased_df_lin_PR <- data.frame(udlPR_goi, udlPR_r2, udlPR_p, stringsAsFactors = FALSE) #create final df
  names(unbiased_df_lin_PR)[1] <- "GOI"
  names(unbiased_df_lin_PR)[2] <- "Rsquared"
  names(unbiased_df_lin_PR)[3] <- "P"
  unbiased_df_lin_PR_sansNA <- na.omit(unbiased_df_lin_PR) #drop anything with R2 or P with NaN
}
unbiased_df_lin_PR$Rsquared[is.nan(unbiased_df_lin_PR$Rsquared)] <- 0
unbiased_df_lin_PR$P[is.nan(unbiased_df_lin_PR$P)] <- 0

#Top
ggplot(unbiased_df_lin_PR, aes(x=Rsquared)) +
  theme_minimal() +
  geom_histogram(binwidth=0.01, fill = "steelblue") +
  xlab("R2") +
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="R2") 

#Bottom
ggplot(unbiased_df_lin_PR, aes(x=P)) +
  theme_minimal() +
  geom_histogram(binwidth=0.01, fill = "steelblue") +
  xlab("P Value") +
  ylab("Frequency") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))  +
  labs(title="P")

###GSEA using Singer 2016 (doi: 10.1016/j.cell.2016.08.052)
T_Modules_Signatures_Raw <- read.csv("/home/a/andramos/t1data/sc_rnaseq/analysis/R_analysis/Analysis/regev2016_genemodule/ref/cpav/T_Modules_Signatures.csv")
#extract equivalet human gene names 
T_Module_Original_NaiveMem <- list('Naive_Mem' = T_Modules_Signatures_Raw$Naïve_Memory_like_module_100)
T_Module_Original_NaiveMem_Features <- na.omit(T_Module_Original_NaiveMem$Naive_Mem)
GenesThatNeedConversionNaiveMem <- T_Module_Original_NaiveMem_Features
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)}
ConvertedGenesNaiveMem <- convertMouseGeneList(GenesThatNeedConversionNaiveMem)

T_Modules_Signatures_Raw <- list('Activation' = T_Modules_Signatures_Raw$Activation_module_edit,
                                 'Activation_Dysfunction' = T_Modules_Signatures_Raw$Activation_Dysfunction_module_edit,
                                 'Extended_Dysfunction' = T_Modules_Signatures_Raw$Dysfunction_module_edit_ext,
                                 'Dysfunction' = T_Modules_Signatures_Raw$Dysfunction_module_edit)
'%!in%' <- function(x,y)!('%in%'(x,y))
T_Modules_Signatures_Dysfunction_Features <- na.omit(T_Modules_Signatures_Raw$Dysfunction)
T_Modules_Signatures_Extended_Dysfunction_Features <- na.omit(T_Modules_Signatures_Raw$Extended_Dysfunction)
T_Modules_Signatures_Activation_Features <- na.omit(T_Modules_Signatures_Raw$Activation)
T_Modules_Signatures_Activation_Dysfunction_Features <- na.omit(T_Modules_Signatures_Raw$Activation_Dysfunction)
T_Modules_Signatures_Naive_Memory_Features <- na.omit(ConvertedGenesNaiveMem)

gene.sets <- list(
  Dysfunction=T_Modules_Signatures_Extended_Dysfunction_Features,
  Activation=T_Modules_Signatures_Activation_Features,
  NaiveMemory=T_Modules_Signatures_Naive_Memory_Features)

ES <- enrichIt(obj = subsetPR, 
               gene.sets = gene.sets, 
               method = "UCell",
               groups = 1000, cores = 4, 
               min.size = NULL)
subsetPR <- AddMetaData(subsetPR, ES)

###Figure4a
##Left
#Top - Dysfunction
VlnPlot(subsetPR, features = c("Dysfunction"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                    "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                    "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                    "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0.01) +  
  NoLegend()  

#Middle - Activation
VlnPlot(subsetPR, features = c("Activation"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                    "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                    "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                    "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0.01) +  
  NoLegend()

#Bottom - Memory
VlnPlot(subsetPR, features = c("NaiveMemory"), 
        same.y.lims = T, pt.size = 0 , cols = c("1-10" =  "#440154FF","11-20" = "#482878FF","21-30" = "#3E4A89FF",
                                                    "31-40" = "#31688EFF","41-50" = "#26828EFF","51-60" = "#1F9E89FF",
                                                    "61-70" = "#35B779FF","71-80" = "#6DCD59FF","81-90" = "#B4DE2CFF",
                                                    "91-100" = "#FDE725FF")) +
  geom_boxplot(color="white", width=0.075, outlier.size = 0.01) +  
  NoLegend()

##Right
UCellEnriched_GSEA <- FetchData(subsetPR, vars =c("Dysfunction","Activation","NaiveMemory", "D52NCAR"))
UCellEnriched_GSEA <- tibble::rownames_to_column(UCellEnriched_GSEA, "barcode")
avg_PR <- function(gene_set){
  gset <- data.frame(matrix(ncol = (length(colnames(gene_set))-1)))
  gset = gset[-1,]
  a1 <- dplyr::semi_join(gene_set, p1, by = "barcode") %>% summarise_if(is.numeric, mean)
  a2 <- dplyr::semi_join(gene_set, p2, by = "barcode") %>% summarise_if(is.numeric, mean)
  a3 <- dplyr::semi_join(gene_set, p3, by = "barcode") %>% summarise_if(is.numeric, mean)
  a4 <- dplyr::semi_join(gene_set, p4, by = "barcode") %>% summarise_if(is.numeric, mean)
  a5 <- dplyr::semi_join(gene_set, p5, by = "barcode") %>% summarise_if(is.numeric, mean)
  a6 <- dplyr::semi_join(gene_set, p6, by = "barcode") %>% summarise_if(is.numeric, mean)
  a7 <- dplyr::semi_join(gene_set, p7, by = "barcode") %>% summarise_if(is.numeric, mean)
  a8 <- dplyr::semi_join(gene_set, p8, by = "barcode") %>% summarise_if(is.numeric, mean)
  a9 <- dplyr::semi_join(gene_set, p9, by = "barcode") %>% summarise_if(is.numeric, mean)
  a10 <- dplyr::semi_join(gene_set, p10, by = "barcode") %>% summarise_if(is.numeric, mean)
  gset <- rbind(gset, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
  gset$Percentile <- c("Percentile 1-10","Percentile 11-20","Percentile 21-30","Percentile 31-40",
                       "Percentile 41-50","Percentile 51-60","Percentile 61-70","Percentile 71-80",
                       "Percentile 81-90","Percentile 91-100")
  gset <- gset %>% dplyr::select(Percentile, everything())
  print("success")
  return(gset)
}
UCEG_PRavg <- avg_PR(UCellEnriched_GSEA)

S16 <- c("Dysfunction","Activation","NaiveMemory")
S16 <- t(UCEG_PRavg[S16])
UBPW_PRN <- function(UBaveraged) {
  resc_df <- NULL
  rowNumbs <- length(rownames(UBaveraged))
  for (i in 1:rowNumbs){
    require(scales)
    resc_row  <- scales::rescale(as.numeric(UBaveraged[i,]))
    resc_df <- rbind(resc_df, resc_row)
  }
  rownames(resc_df) <- rownames(UBaveraged)
  colnames(resc_df) <- colnames(UBaveraged)
  resc_df <- as.data.frame(resc_df, stringsAsFactors = FALSE)
  print("success")
  return(resc_df)
}
S16n <- UBPW_PRN(S16)
S16n <- t(S16n)
S16n <- data.frame(S16n, UCEG_PRavg$D52NCAR)
colnames(S16n)[4] <- "D52NCAR"
rownames(S16n) <- c("1-10", "11-20","21-30", "31-40","41-50",
                    "51-60","61-70","71-80","81-90","91-100")
#Top - Dysfunction
S16d_log <- lm(Dysfunction ~ log(D52NCAR), data = S16n)
ggplot(S16n, mapping = aes(x=D52NCAR, y=Dysfunction)) +
  geom_point(size =2) + theme_linedraw() +
  theme(text=element_text(size=15)) +
  stat_smooth(method='lm', formula = y ~ log(x), se=F) +
  geom_text_repel(label=UCEG_PRavg$Percentile, size = 6) +
  labs(x="",y = "",
       subtitle = paste("Radj2=",signif(summary(S16d_log)$adj.r.squared,3),
                        ", P =", signif(summary(S16d_log)$coefficients[2,4],3))) +
  scale_y_continuous(limits=c(-0.135, 1.125), (breaks=seq(0,1.125,0.125)))

#Middle - Activation
S16a_log <- lm(Activation ~ log(D52NCAR), data = S16n)
ggplot(S16n, mapping = aes(x=D52NCAR, y=Activation)) +
  geom_point(size = 2) + theme_linedraw() +
  theme(text=element_text(size=15)) +
  stat_smooth(method='lm', formula = y ~ log(x), se = F) +
  geom_text_repel(label=UCEG_PRavg$Percentile, size = 6) +
  labs(x="",y = "",
       subtitle = paste("Radj2=",signif(summary(S16a_log)$adj.r.squared,3),
                        ", P =", signif(summary(S16a_log)$coefficients[2,4],3))) +
  scale_y_continuous(limits=c(-0.135, 1.125), (breaks=seq(0,1.125,0.125)))

#Bottom - Memory
S16nm_log <- lm(NaiveMemory ~ log(D52NCAR), data = S16n)
ggplot(S16n, mapping = aes(x=D52NCAR, y=NaiveMemory)) +
  geom_point(size = 2) + theme_linedraw() +
  theme(text=element_text(size=15)) +
  stat_smooth(method='lm', formula = y ~ log(x), se = F) +
  geom_text_repel(label=UCEG_PRavg$Percentile, size = 6) +
  labs(x="",y = "",
       subtitle = paste("Radj2=",signif(summary(S16nm_log)$adj.r.squared,3),
                        ", P =", signif(summary(S16nm_log)$coefficients[2,4],3))) +
  scale_y_continuous(limits=c(-0.135, 1.125), (breaks=seq(0,1.125,0.125)))

###Extended Data Figure 15
T_Modules_Signatures_Extended_Dysfunction_Features_CAR <- append("D52NCAR", T_Modules_Signatures_Extended_Dysfunction_Features)
T_Modules_Signatures_Activation_Features_CAR <- append("D52NCAR", T_Modules_Signatures_Activation_Features)
T_Modules_Signatures_Naive_Memory_Features_CAR <- append("D52NCAR", T_Modules_Signatures_Naive_Memory_Features)
DysGenes <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% T_Modules_Signatures_Extended_Dysfunction_Features_CAR,], stringsAsFactors = F)
ActGenes <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% T_Modules_Signatures_Activation_Features_CAR,], stringsAsFactors = F)
NMemGenes <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% T_Modules_Signatures_Naive_Memory_Features_CAR,], stringsAsFactors = F)

#Dysfunction - Top
rownames(DysGenes) <- DysGenes$GOI
DysGenes <- DysGenes[,-1]
DysGenes <- na.omit(DysGenes)
colnames(DysGenes) <- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(DysGenes), cluster_cols = F)
#Activation - Middle
rownames(ActGenes) <- ActGenes$GOI
ActGenes <- ActGenes[,-1]
colnames(ActGenes) <- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(ActGenes), cluster_cols = F)
#Memory - Right
rownames(NMemGenes) <- NMemGenes$GOI
NMemGenes <- NMemGenes[,-1]
colnames(NMemGenes)<- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(NMemGenes), cluster_cols = F)

###Unbiased exploration analysis
##Mutual Information
#Figure 4b

UBPRMGenesa_T <- as.data.frame(t(UBPRMGenesa), stringsAsFactors = FALSE)
  
CARdisc <- discretize(UBPRMGenesa_T$D52NCAR)
CARE <- infotheo::entropy(CARdisc)
calc_MI <- NULL
  for (i in 1:(which(colnames(UBPRMGenesa_T) == "D52NCAR"))) {
    GOIdisc <- discretize(UBPRMGenesa_T[,i])
    GOIE <- infotheo::entropy(GOIdisc)
    MI <- mutinformation(CARdisc, GOIdisc)
    NMI <- (MI / sqrt(CARE*GOIE))
    a <- cbind.data.frame(MI,NMI)
    calc_MI <- rbind.data.frame(calc_MI, a)
  }
GOI_NMI <- calc_NMI
rownames(GOI_NMI) <- colnames(UBPRMGenesa_T)
GOI_NMI <- na.omit(GOI_NMI)
GOI_MaxMI <- dplyr::filter(GOI_NMI, MI == max((GOI_NMI$MI)))

MI1GOI <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% rownames(GOI_MaxMI),], stringsAsFactors = F)
rownames(MI1GOI) <- MI1GOI$GOI
MI1GOI <- MI1GOI[,-1]
colnames(MI1GOI) <- c("1-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(MI1GOI), cluster_cols = F, clustering_method = "mcquitty")

#Extended Data Figure 16 
#Linear regression 
car_index <- which(rownames(UBPRMGenesa) == "D52NCAR")  #CAR is the last row
length_sc <- length(rownames(UBPRMGenesa))-1 #subtracted one as our CAR is the last row
unbiased_df_lin_PR <- NULL                 #initialise empty df to fill                                       
  for (i in 1:length_sc) {
    ub_linfit_PR <- lm(unlist(UBPRMGenesa[i,]) ~ unlist(UBPRMGenesa[15389,])) #conduct linear fit
    goi_PR <- rownames(UBPRMGenesa[i,]) #gene name
    r2_PR <- as.numeric(summary(ub_linfit_PR)$r.squared) #Rsquared
    pval_PR <- as.numeric(summary(ub_linfit_PR)$coefficients[2,4]) #p value
    linf_PR <- cbind(goi_PR, r2_PR, pval_PR) #create new row
    unbiased_df_lin_PR <- rbind(unbiased_df_lin_PR, linf_PR) #add to initialised df 
  }                                 
unbiased_df_lin_PR <- as.data.frame(unbiased_df_lin_PR, stringsAsFactors=FALSE) #create intermediate df
udlPR_goi <- unbiased_df_lin_PR$goi
udlPR_r2 <- unbiased_df_lin_PR$r2 #convert from character to numeric
udlPR_r2 <- as.numeric(udlPR_r2)
udlPR_p <- unbiased_df_lin_PR$pval
udlPR_p <- as.numeric(udlPR_p) #convert from character to numeric
unbiased_df_lin_PR <- data.frame(udlPR_goi, udlPR_r2, udlPR_p, stringsAsFactors = FALSE) #create final df
names(unbiased_df_lin_PR)[1] <- "GOI"
names(unbiased_df_lin_PR)[2] <- "Rsquared"
names(unbiased_df_lin_PR)[3] <- "P"
unbiased_df_lin_PR_sansNA <- na.omit(unbiased_df_lin_PR) #drop anything with R2 or P with NaN

UB_LinearRegression <- unbiased_df_lin_PR_sansNA %>% arrange(desc(Rsquared)) 

#Logarithmic regression 
car_index <- which(rownames(UBPRMGenesa) == "D52NCAR") 
length_sc <- length(rownames(UBPRMGenesa))-1
unbiased_df_log_PR <- NULL                                                       
  for (i in 1:length_sc) {
    ub_logfit_PR <- lm(unlist(UBPRMGenesa[i,]) ~ unlist(log(UBPRMGenesa[15389,])))
    goi_PR <- rownames(UBPRMGenesa[i,])
    r2_PR <- as.numeric(summary(ub_logfit_PR)$r.squared)
    pval_PR <- as.numeric(summary(ub_logfit_PR)$coefficients[2,4])
    logf_PR <- cbind(goi_PR, r2_PR, pval_PR)
    unbiased_df_log_PR <- rbind(unbiased_df_log_PR, logf_PR)
  }                                 
unbiased_df_log_PR <- as.data.frame(unbiased_df_log_PR, stringsAsFactors=FALSE)
udloPR_goi <- unbiased_df_log_PR$goi
udloPR_r2 <- unbiased_df_log_PR$r2
udloPR_r2 <- as.numeric(udloPR_r2)
udloPR_p <- unbiased_df_log_PR$pval
udloPR_p <- as.numeric(udloPR_p)
unbiased_df_log_PR <- data.frame(udloPR_goi, udloPR_r2, udloPR_p, stringsAsFactors = FALSE)
names(unbiased_df_log_PR)[1] <- "GOI"
names(unbiased_df_log_PR)[2] <- "Rsquared"
names(unbiased_df_log_PR)[3] <- "P"
unbiased_df_log_PR_sansNA <- na.omit(unbiased_df_log_PR)

UB_LogarithmicRegression <- unbiased_df_log_PR_sansNA %>% arrange(desc(Rsquared)) 

#EDF16a
Metab_Genes <- read_excel("/Users/ahsramos/Desktop/Manuscript/2023_Drafting/Figures/Metabolism/GO.xlsx")
KOP <- Metab_Genes$OxPhos
KOP <- append(KOP, "D52NCAR")
KOxPhos <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% KOP,], stringsAsFactors = F)
rownames(KOxPhos) <- KOxPhos$GOI
KOxPhos <- KOxPhos[,-1]
colnames(KOxPhos) <- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(KOxPhos), cluster_cols = F)

#EDF16e
#GO:0061621 - Amigo Canonical Glycolysis
canonglyc <- c("PKM","TPI1","GCK","PFKP", "PGAM1","PGAM2", "PFKM", "FOXK1", "ENO1",
               "ENO3","HK3", "HK2", "FOXK2", "PFKL", "ENO2", "HK1", "PGK1", "PGK2",
               "HEL-S-30","HEL-S-68p","HEL-S-49")
CanonGly <- as.data.frame(UBPRMGenesa_zed[UBPRMGenesa_zed$GOI %in% canonglyc ,], stringsAsFactors = F)
rownames(CanonGly) <- CanonGly$GOI
CanonGly <- CanonGly[,-1]
colnames(CanonGly) <- c("1-10", "11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100")
pheatmap(as.matrix(na.omit(CanonGly)), cluster_cols = F)
