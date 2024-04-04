library(ATACseqQC)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(ChIPpeakAnno)
library(Rsamtools)

bamFile="wt_a4_4h_1.rmChrM.bam"
bamFileLabels <- "wt_a4_4h_1"

fragSize <- fragSizeDist(bamFile, bamFileLabels)

bam_qc=bamQC(bamFile, outPath = NULL)
bam_qc[1:10]
outPath <- "splitBam"
dir.create(outPath)

seqlev <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
which <- as(seqinfo(Athaliana)[seqlev], "GRanges")

gal <- readBamFile(bamFile, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal)

txs <- transcripts(TxDb.Athaliana.BioMart.plantsmart28)
newStyle <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC")
txs <- renameSeqlevels(txs, newStyle)

pt <- PTscore(gal1, txs)
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")

nfr <- NFRscore(gal1, txs)
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))

tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore

plot(100*(-9:10-.5), tsse$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")

genome <- Athaliana
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)
dir(outPath)

bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
## estimate the library size for normalization
(librarySize <- estLibSize(bamfiles))

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  ylab="Averaged coverage")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")



QC_function <- function(bamFile, bamFileLabels) {
  fragSize <- fragSizeDist(bamFile, bamFileLabels)
  
  bam_qc=bamQC(bamFile, outPath = NULL)
  print(bam_qc[1:10])

  seqlev <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
  which <- as(seqinfo(Athaliana)[seqlev], "GRanges")

  gal <- readBamFile(bamFile, which=which, asMates=TRUE, bigFile=TRUE)
  gal1 <- shiftGAlignmentsList(gal)

  txs <- transcripts(TxDb.Athaliana.BioMart.plantsmart28)
  newStyle <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrM", "ChrC")
  txs <- renameSeqlevels(txs, newStyle)

  pt <- PTscore(gal1, txs)
  ptplot <- plot(pt$log2meanCoverage, pt$PT_score,
                 xlab="log2 mean coverage",
                 ylab="Promoter vs Transcript")
  print(ptplot)

  nfr <- NFRscore(gal1, txs)
  nfrplot <- plot(nfr$log2meanCoverage, nfr$NFR_score,
                  xlab="log2 mean coverage",
                  ylab="Nucleosome Free Regions score",
                  main="NFRscore for 200bp flanking TSSs",
                  xlim=c(-10, 0), ylim=c(-5, 5))
  print(nfrplot)

  tsse <- TSSEscore(gal1, txs)
  tsse$TSSEscore

  tssplot <- plot(100*(-9:10-.5), tsse$values, type="b",
                  xlab="distance to TSS",
                  ylab="aggregate TSS score")
  print(tssplot)

  genome <- Athaliana
  objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = "splitBam")

  bamfiles <- file.path("splitBam",
                        c("NucleosomeFree.bam",
                          "mononucleosome.bam",
                          "dinucleosome.bam",
                          "trinucleosome.bam"))
  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  ## estimate the library size for normalization
  (librarySize <- estLibSize(bamfiles))

  NTILE <- 101
  dws <- ups <- 1010
  sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                       "mononucleosome",
                                       "dinucleosome",
                                       "trinucleosome")],
                            TSS=TSS,
                            librarySize=librarySize,
                            seqlev=seqlev,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)
  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
  fAheatmap <- featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                                     zeroAt=.5, n.tile=NTILE)
  print(fAheatmap)

  out <- featureAlignedDistribution(sigs,
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, n.tile=NTILE, type="l",
                                    ylab="Averaged coverage")

  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  out <- apply(out, 2, range01)
  posplot <- matplot(out, type="l", xaxt="n",
                     xlab="Position (bp)",
                     ylab="Fraction of signal")
  axis(1, at=seq(0, 100, by=10)+1,
       labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
  abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
  print(posplot)
}

bams <- list.files(pattern="\\.bam$")
bamLabels <- list.files(pattern="\\.bam.bai$")

for (i in c(1,2)) {
  QC_function(bams[i], bamLabels[i]) 
}
