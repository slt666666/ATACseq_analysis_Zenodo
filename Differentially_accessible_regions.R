library(edgeR)
library(EDASeq)

library(GenomicAlignments)
library(GenomicFeatures)

library(dplyr)

ff=FaFile("TAIR10_chr_all.fas")

treats_set=list(c("wt_a4_4h", "wt_mk_4h"), c("wt_kv_4h", "wt_mk_4h"), c("setiwt_est", "setiwt_mk"), c("setikv_est", "setikv_mk"))
for (treats in treats_set) {
  # Arrange data table format
  cnt_table = read.table(paste0(treats[1], "_", treats[2], "_merged_peaks.counts.tsv"), sep="\t", header=TRUE, blank.lines.skip=TRUE)
  rownames(cnt_table)=cnt_table$Geneid
  sample_names = colnames(cnt_table)[7:dim(cnt_table)[2]]
  new_colnames=c("Geneid","Chr","Start","End","Strand","Length")
  groups = c()
  for (sample_name in sample_names){
    new_colnames=c(new_colnames, substr(sample_name, 17, nchar(sample_name)-17))
    groups=c(groups, substr(sample_name, 17, nchar(sample_name)-19))
  }
  colnames(cnt_table)=new_colnames
  reads.peak = cnt_table[,c(7:dim(cnt_table)[2])]
  
  # Prepare GC content data for GC-aware normalization
  gr = GRanges(seqnames=cnt_table$Chr, ranges=IRanges(cnt_table$Start, cnt_table$End), strand="*", mcols=data.frame(peakID=cnt_table$Geneid))
  peakSeqs = getSeq(x=ff, gr)
  gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
  
  # calculate the offsets which correct for library size as well as GC content
  reads.peak=as.matrix(reads.peak)
  dataOffset = withinLaneNormalization(reads.peak,y=gcContentPeaks,num.bins=20,which="full",offset=TRUE)
  dataOffset = betweenLaneNormalization(reads.peak,which="full",offset=TRUE)
  
  # Identify differential accessibility
  design = model.matrix(~groups)
  
  d = DGEList(counts=reads.peak, group=groups)
  
  keep = filterByExpr(d)
  d=d[keep,,keep.lib.sizes=FALSE]
  
  d$offset = -dataOffset[keep,]
  d.eda = estimateGLMCommonDisp(d, design = design)
  d.eda = estimateGLMCommonDisp(d, design = design)
  fit = glmFit(d.eda, design = design)
  lrt.EDASeq = glmLRT(fit, coef = 2)
  
  DA_res=as.data.frame(topTags(lrt.EDASeq, nrow(lrt.EDASeq$table)))
  DA_res$Geneid = rownames(DA_res)
  DA.res.coords = left_join(DA_res,cnt_table[1:4],by="Geneid")
  
  write.table(DA.res.coords, paste0(treats[1], "_vs_", treats[2], ".tsv"), quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, fileEncoding = "")
}