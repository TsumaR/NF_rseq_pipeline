#-- process start --#

#-- read cleaning by flexbar --#
$flexbar -z GZ -n $cpu -u 10 -a $adapter -at ANY --htrim-right AT \
  --htrim-min-length 10 --htrim-error-rate 0.1 -t $cln_pref \
  -r $raw1 -p $raw2
#-- read cleaning by fastq-mcf(ea-utils) --#
#  $fastq_mcf \
#    -l 25 -x 10 \
#    --qual-mean 20 \
#    -H --homopolymer-pct 50 \
#    -o $cln1 -o $cln2 \
#    $adapter $raw1 $raw2
#  gzip $cln1
#  gzip $cln2

#-- align to genome --#
$hisat2 -p $cpu -x $ht2_index -1 $cln1gz -2 $cln2gz | \
  $samtools view -b | $samtools sort -@ $cpu > $aln

#-- region tagging --#
$samtools view -b -f4 $aln > $aln_unmap
$samtools view -h -F4 $aln | grep -v $'\t'ERCC | $samtools view -b > $aln_hits
$samtools view -h -F4 $aln | grep -e '^@' -e $'\t'ERCC | $samtools view -b > $aln_ercc
$bedtools intersect -wa -u -a $aln_hits -b $bed_gene > $aln_gene
$bedtools intersect -wa -v -a $aln_hits -b $bed_gene > $aln_intergenic
$bedtools intersect -wa -u -a $aln_gene -b $bed_exon > $aln_exon
$bedtools intersect -wa -v -a $aln_gene -b $bed_exon > $aln_intron
$region_tagger -o $aln_exon $aln_exon exon
$region_tagger -o $aln_intron $aln_intron intron
$region_tagger -o $aln_intergenic $aln_intergenic intergenic
$region_tagger -o $aln_ercc $aln_ercc ercc
$samtools merge -f -@ $cpu $aln $aln_exon $aln_intron $aln_intergenic $aln_ercc
$samtools sort -@ $cpu -o $aln $aln
$samtools index $aln

#-- estimate gene expression --#
$stringtie -e -p $cpu -A $gene_abund -o $trans_gtf -G $gtf $aln

#-- isoform expression --#
$extract_isoform_expression $trans_gtf > $iso_abund

# coverage
$bam2wig -i $aln -s $genome_size -o $aln_pref

#-- quality check --#
$fastqc $raw1
$fastqc $raw2
$fastqc $cln1gz
$fastqc $cln2gz

#-- qc for R --#
$qc_alignment -o $qc_out_aln1 $fqc_raw1 $fqc_cln1 $aln
$qc_alignment --read2 -o $qc_out_aln2 $fqc_raw2 $fqc_cln2 $aln
$qc_distribution -o $qc_out_dist1 $fqc_cln1 $aln
$qc_distribution --read2 -o $qc_out_dist2 $fqc_cln2 $aln
$qc_insert_size -o $qc_out_isize $aln
$qc_coverage --isoform $iso_abund -o $qc_out_cov $gtf ${aln_pref}.bw

#-- remove intermediate files --#
rm ${aln_pref}.*.bam ${aln_pref}.wig

