#-- process start --#

#-- read cleaning --#
$flexbar -z GZ -n $cpu -u 10 -a $adapter -at ANY --htrim-right AT \
  --htrim-min-length 10 --htrim-error-rate 0.1 -t $cln_pref \
  -r $raw

#-- align to genome --#
$hisat2 -p $cpu -x $ht2_index -U $clngz | \
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
$fastqc $raw
$fastqc $clngz

#-- qc for R --#
$qc_alignment -o $qc_out_aln $fqc_raw $fqc_cln $aln
$qc_distribution -o $qc_out_dist $fqc_cln $aln
$qc_coverage --isoform $iso_abund -o $qc_out_cov $gtf ${aln_pref}.bw

#-- remove intermediate files --#
rm ${aln_pref}.*.bam ${aln_pref}.wig

