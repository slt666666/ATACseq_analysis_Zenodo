name=sample_name

# mapping by bowtie2
bowtie2 --threads 8 \
          -x ../reference2/bowtie2_index \
          --very-sensitive \
          --end-to-end \
          --maxins 20000 \
          --no-discordant \
          --no-mixed \
          --time \
          --no-unal \
          --qc-filter \
          -1 ${name}_R1.fastq.gz \
          -2 ${name}_R2.fastq.gz | \
  samtools view -b --reference TAIR10_chr_all.fas --threads 8 | \
  samtools sort -o ${name}.sort.bam && samtools index ${name}.sort.bam

# remove reads mapped to mitochondria & chloroplast genome
samtools view -h ${name}.sort.bam | grep -v ChrM | grep -v ChrC | samtools sort -O bam -o ${name}.rmChrM.bam -T .
samtools index ${name}.rmChrM.bam

# Mark duplicates & Remove duplicates & low-quality alignments
gatk MarkDuplicates QUIET=true INPUT=${name}.rmChrM.bam OUTPUT=${name}.marked.bam METRICS_FILE=${name}.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
samtools view -h -b -f 2 -F 1548 -q 30 ${name}.marked.bam | samtools sort -o ${name}.filtered.bam
samtools index ${name}.filtered.bam

# Shift read coordinates
alignmentSieve --numberOfProcessors max --ATACshift --bam ${name}.filtered.bam -o ${name}.shifted.bam
samtools sort -O bam -o ${name}.shifted.sort.bam ${name}.shifted.bam
samtools index ${name}.shifted.sort.bam

# Peak calling by MACS3
macs3 callpeak -f BAMPE --call-summits -c gdna.shifted.bam -t ${name}.shifted.sort.bam -g 119482012 -n ${name}.macs3.summits.bampe -B -q 0.05

# Get consensus peaks for each treatment
for treat in wt_a4_4h wt_kv_4h wt_mk_4h wt_un setikv_est setikv_mk setikv_un setiwt_est setiwt_mk setiwt_un; do
  array=($(ls ${treat}*.macs3.summits.bampe_peaks.narrowPeak))
  if [ 2 -eq ${#array[@]} ]; then
    bedtools intersect -a ${array[0]} -b ${array[1]} -f 0.50 -r > tmp.bed
  else
    bedtools intersect -a ${array[0]} -b ${array[1]} ${array[2]} -f 0.50 -r > tmp.bed
  fi
  awk '! a[$1" "$2" "$3" "$4]++' tmp.bed > ${treat}_peaks.bed
done

# Count the number of reads in peaks for the investigation of differential accesiblity
treat1=wt_a4_4h
treat2=wt_mk_4h
bedops -m ${treat1}_peaks.bed ${treat2}_peaks.bed > ${treat1}_${treat2}_merged_peaks.bed
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="wt_a4_mk_merge_"++nr;  print peakid,$1,$2,$3,"."}' ${treat1}_${treat2}_merged_peaks.bed > ${treat1}_${treat2}_merged_peaks.saf
featureCounts -p -F SAF -a ${treat1}_${treat2}_merged_peaks.saf --fracOverlap 0.2 \
              -o ${treat1}_${treat2}_merged_peaks.counts \
              ${treat1}*shifted.sort.bam ${treat2}*shifted.sort.bam