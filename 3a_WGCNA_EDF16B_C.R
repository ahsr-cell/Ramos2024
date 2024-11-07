### WGCNA
require(tidyverse)     
require(magrittr)      
require(WGCNA)
require(rbioapi)
require(gghighlight)
require(ggVennDiagram)
setwd("/Users/ahsramos/Desktop/WGCNA")

##selecting threshold values
allowWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  UBPRMGenesa_T,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#based on these visualisations, choose 7-9 as we want a value close to the curve 

picked_power = 8
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)

#Conduct WGCNA for high resolution - use mergeCutHeight = 0.25
netwk_high_resolution <- blockwiseModules(UBPRMGenesa_T,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          #for pure positive and negative correlation, used 0.25 
                          #for low resolution, used 0.45
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

#Conduct WGCNA for low resolution - use mergeCutHeight = 0.45
netwk_low_resolution <- blockwiseModules(UBPRMGenesa_T,                # <= input here
                                          
                                          # == Adjacency Function ==
                                          power = picked_power,                # <= power here
                                          networkType = "signed",
                                          
                                          # == Tree and Block Options ==
                                          deepSplit = 2,
                                          pamRespectsDendro = F,
                                          # detectCutHeight = 0.75,
                                          minModuleSize = 30,
                                          maxBlockSize = 4000,
                                          
                                          # == Module Adjustments ==
                                          reassignThreshold = 0,
                                          mergeCutHeight = 0.45,
                                          #for pure positive and negative correlation, used 0.25 
                                          #for low resolution, used 0.45
                                          # == TOM == Archive the run results in TOM file (saves time)
                                          saveTOMs = T,
                                          saveTOMFileBase = "ER",
                                          
                                          # == Output Options
                                          numericLabels = T,
                                          verbose = 3)
cor <- temp_cor
## Convert labels to colors for plotting
#Highresolution
mergedColors_high_resolution = labels2colors(netwk_high_resolution$colors)
#Low resolution
mergedColors_low_resolution = labels2colors(netwk_low_resolution$colors)

## Plot the dendrogram and the module colors underneath
#High resolution
plotDendroAndColors(
  netwk_high_resolution$dendrograms[[1]],
  mergedColors[netwk_high_resolution$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#Low resolution
plotDendroAndColors(
  netwk_low_resolution$dendrograms[[1]],
  mergedColors[netwk_low_resolution$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

#High resolution
module_df_high_resolution <- data.frame(
  gene_id = names(netwk_high_resolution$colors),
  colors = labels2colors(netwk_high_resolution$colors)
)

#Low resolution
module_df_low_resolution <- data.frame(
  gene_id = names(netwk_low_resolution$colors),
  colors = labels2colors(netwk_low_resolution$colors)
)

# Get Module Eigengenes per cluster
#High resolution
MEs0_high_resolution <- moduleEigengenes(UBPRMGenesa_T, mergedColors_high_resolution)$eigengenes
#Low resolution
MEs0_low_resolution <- moduleEigengenes(UBPRMGenesa_T, mergedColors_low_resolution)$eigengenes

# Reorder modules so similar modules are next to each other
#High resolution
MEs0_high_resolution <- orderMEs(MEs0_high_resolution)
module_order_high_resolution = names(MEs0_high_resolution) %>% gsub("ME","", .)

#Low resolution
MEs0_low_resolution <- orderMEs(MEs0_low_resolution)
module_order_low_resolution = names(MEs0_low_resolution) %>% gsub("ME","", .)

# Add treatment names
#High resolution
MEs0_high_resolution$treatment = row.names(MEs0_high_resolution)

#Low resolution
MEs0_low_resolution$treatment = row.names(MEs0_low_resolution)

# tidy & plot data
mME_high_resolution = MEs0_high_resolution %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order_high_resolution)
  )

mME_low_resolution = MEs0_low_resolution %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order_low_resolution )
  )

###Extended Data Figure 16c, top-left
mME_high_resolution %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "High Resolution: Module-CAR Expression Relationships",x = "CAR Expression Bin", y = "Modules", fill="corr")

###Extended Data Figure 16c, bottom-left 
mME_low_resolution %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Low Resolution: Module-CAR Expression Relationships",x = "CAR Expression Bin", y = "Modules", fill="corr")

#Modules of interest
modules_of_interest_high_resolution = c("turquoise", "brown")
modules_of_interest_low_resolution = c("magenta")

brown_mod <- filter(module_df_high_resolution, colors == "brown")
turquoise_mod <- filter(module_df_high_resolution, colors == "turquoise")
magenta_mod <- filter(module_df_low_resolution, colors == "magenta")

#write_csv(brown_mod,"brown.csv")
#write_csv(turquoise_mod,"turquoise.csv")
#write_csv(magenta_mod,"magenta.csv")

# Pull out list of genes in that module

modcount <- filter(module_df, colors == "turquoise")
length(modcount$gene_id)

###Extended Data Figure 16c, top-right 
modules_of_interest = c("turquoise")
submod = module_df_high_resolution %>%
  subset(colors %in% modules_of_interest)

row.names(module_df_high_resolution) = module_df_high_resolution$gene_id

UBPRMGenesa_zed <- column_to_rownames(UBPRMGenesa_zed, "GOI") #execute if GOI is its own column still 
subexpr =UBPRMGenesa_zed[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_high_resolution[gene_id,]$colors
  )
Turq_Average <- colMeans(subexpr)
Turq_Average <- as.data.frame(Turq_Average)
names(Turq_Average)[1] <- "Average"
Turq_Average <- t(Turq_Average)
subexpr_turquoise <- subexpr
subexpr_turquoise <- rbind(subexpr_turquoise, Turq_Average)

submod_df_turq = data.frame(subexpr_turquoise) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_high_resolution[gene_id,]$colors
  )

submod_df_turq %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "turquoise"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score")

#Isolate Average line
submod_df_turq %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "turquoise"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score") + 
  gghighlight( gene_id == "Average"
              ) 

###Extended Data Figure 16c, middle-right
modules_of_interest = c("magenta")
submod = module_df_low_resolution %>%
  subset(colors %in% modules_of_interest)

row.names(module_df_low_resolution) = module_df_low_resolution$gene_id

subexpr =UBPRMGenesa_zed[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_low_resolution[gene_id,]$colors
  )

Mag_Average <- colMeans(subexpr)
Mag_Average <- as.data.frame(Mag_Average)
names(Mag_Average)[1] <- "Average"
Mag_Average <- t(Mag_Average)
subexpr_magenta <- subexpr
subexpr_magenta <- rbind(subexpr_magenta, Mag_Average)

submod_df_mag = data.frame(subexpr_magenta) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_low_resolution[gene_id,]$colors
  )

submod_df_mag %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "magenta"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score")

#Isolate Average line
submod_df_mag %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "magenta"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score") + 
  gghighlight( gene_id == "Average"
  ) 

###Extended Data Figure 16c, bottom-right 
modules_of_interest = c("brown")
submod = module_df_high_resolution %>%
  subset(colors %in% modules_of_interest)

subexpr =UBPRMGenesa_zed[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_high_resolution[gene_id,]$colors
  )

Brown_Average <- colMeans(subexpr)
Brown_Average <- as.data.frame(Brown_Average)
names(Brown_Average)[1] <- "Average"
Brown_Average <- t(Brown_Average)
subexpr_brown <- subexpr
subexpr_brown <- rbind(subexpr_brown, Brown_Average)

submod_df_br = data.frame(subexpr_brown) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df_high_resolution[gene_id,]$colors
  )

submod_df_br %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "brown"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score") 

#Isolate average line
submod_df_br %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = "brown"),
            alpha = 0.2,  show.legend = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90))  +
  labs(x = "CAR Expression Bin",
       y = "Z-score") + 
  gghighlight( gene_id == "Average"
  ) 


### Extended Data Figure 16b
setwd("/Users/ahsramos/Desktop/WGCNA")
magenta_mod <- read_csv("/Users/ahsramos/Desktop/WGCNA/magenta.csv")
brown_mod <- read_csv("/Users/ahsramos/Desktop//WGCNA/brown.csv")
turquoise_mod <- read_csv("/Users/ahsramos/Desktop/WGCNA/turquoise.csv")

WGCNAgenes <- unique(c(magenta_mod$gene_id, brown_mod$gene_id, turquoise_mod$gene_id))
length(WGCNAgenes) #reports 5066

GOI_MaxMI #check exists - see 2_GeneExpressionMatrix_UnbiasedDiscovery.R, line 412
length(rownames(GOI_MaxMI)) #reports 977, includes D52NCAR, therefore subtract one
UBLin_GOI_strong  <- as.data.frame(filter(UB_LinearRegression, Rsquared >= 0.9)) #pulled from 2_GeneExpressionMatrix_UnbiasedDiscovery.R, line 445
length(UBLin_GOI_strong$GOI) #reports 221, does not include D52NCAR
UBLog_GOI_strong <- as.data.frame(filter(UB_LogarithmicRegression, Rsquared >= 0.9)) #pulled from 2_GeneExpressionMatrix_UnbiasedDiscovery.R, line 471
length(UBLog_GOI_strong$GOI) #reports 365, does not include D52NCAR

ggVennDiagram(VD_NMI_UB_WGCNA, label_alpha = 0, category.names = c("Linear: 221","Logistic: 365","Mutual Information: 977", "WGCNA")) +
  ggplot2::scale_fill_gradient(low="white",high = "red")
