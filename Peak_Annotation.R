library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(clusterProfiler)
library(bedr)
library(UpSetR)
library(diffloop)

input_file_1 <- "wt_a4_4h_1.macs3.summits.bampe_peaks.narrowPeak"
input_file_2 <- "wt_a4_4h_2.macs3.summits.bampe_peaks.narrowPeak"

files <- c(input_file_1, input_file_2)

# make promoter regions & fix seqlevels
txdb <- TxDb.Athaliana.BioMart.plantsmart28
promoter <- getPromoters(TxDb=txdb,upstream=2000,downstream=2000)
seqlevels(promoter) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "mitochondria", "chloroplast")

# Profile of ATAC peaks binding to TSS regions for all treat
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
names(tagMatrixList) <- c("rep1", "rep2")
plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), origin_label="TSS")

# Profile of ATAC peaks binding to TTS regions for all treat
TTS_region <- read.table("TTS.bed")
colnames(TTS_region) <- c("chr", "start", "end", "strand")
TTS_region <- makeGRangesFromDataFrame(TTS_region)
seqlevels(TTS_region) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "mitochondria", "chloroplast")

tagMatrixList_TTS <- lapply(files, getTagMatrix, windows=TTS_region)
names(tagMatrixList_TTS) <- c("rep1", "rep2")
plotAvgProf(tagMatrixList_TTS, xlim=c(-2000, 2000), origin_label="TTS")

# Peak heatmap of ATAC peaks binding to TTS regions for all treat
tagHeatmap(tagMatrixList, xlim=c(-2000, 2000), color=c("#D55E00", "#0072B2", "#E69F00", "black"))

# Visualize Genomic Annotation
peak1 <- readPeakFile(files[1])
seqlevels(peak1) <- c("1", "2", "3", "4", "5")
plotAnnoBar(annotatePeak(peak1, tssRegion=c(-2000, 1000),TxDb=txdb))
peak2 <- readPeakFile(files[2])
seqlevels(peak2) <- c("1", "2", "3", "4", "5")
plotAnnoBar(annotatePeak(peak2, tssRegion=c(-2000, 1000),TxDb=txdb))

# Overlap of peaks and annotated genes
peakAnnoList <- lapply(list(peak1, peak2), annotatePeak, TxDb=txdb,
                       tssRegion=c(-2000, 1000), verbose=FALSE)
genes=lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) <- c("rep1", "rep2")
vennplot(genes)

# visualize by UpSet plot
upset(fromList(genes), order.by = "freq")