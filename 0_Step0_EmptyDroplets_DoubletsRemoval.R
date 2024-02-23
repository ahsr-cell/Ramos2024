### Empty droplet and Doublet/multiplet filtering  filtering 

require(DropletUtils) 
require(scDblFinder)

#EmptyDrops filtering
raw_path <- "/home/a/andramos/t1data/sc_rnaseq/analysis/cellranger/2feb/D52N_miSFITs/outs/raw_feature_bc_matrix"
sce <- read10xCounts(raw_path, col.names=TRUE)

bcrank <- barcodeRanks(sce, lower = 500)
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq],
     bcrank$total[uniq],
     log="xy",
     xlab="Rank",
     ylab="Total UMI count",
     cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

set.seed(500) 

#defining cells with EmptyDrops conditions
e.out <- emptyDrops(counts(sce), lower = 500)
e.keep<- e.out$FDR <= 0.01
summary(e.keep)

#and with Cellranger
c.keep <- defaultDrops(counts(sce))
summary(c.keep)

#retaining cells detected by either EmptyDrops or Cellranger
keep <- c.keep | (e.keep & !is.na(e.keep))
summary(keep)

#Adding to metadata if needed for future use
sce$PValue <- e.out$PValue

#updating barcodes matrix with EmptyDrops results 
##note barcode list varies due to empricism/'randomness' (e.g., setseed) of EmptyDrops algorithm (e.g., #BCs = 20339, 20340, 20350, etc. )
tenx.counts.filtered <- sce[,(keep)]
EDfilteredBarcodes <- data.frame(sort(colnames(tenx.counts.filtered))) #used for qualititative check 

#scDblfinder for in silico doublet/multiplet annotation
##note annotation can slightly vary due to model building/'randomness'
scDbl <- scDblFinder(tenx.counts.filtered)
scDblBarcodes <- data.frame(scDbl@assays@data@listData[["counts"]]@Dimnames[[2]],scDbl$scDblFinder.class) 

#Pure barcode list provided to genetic demultiplexing pipeline - annotations later used in Demultiplexing Converging 
write_tsv(x = scDblBarcodes,
          file = str_glue('/home/a/andramos/ramoscandido/6_tenxpreprocessing/10X-demux/gendem/gPlexA/results.gPlexA.dir/adapting/scDblFinder/scDblBarcodes.tsv'),
          col_names = FALSE)  

