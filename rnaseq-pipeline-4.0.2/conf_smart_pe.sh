#-- references --#
idx_dir=$HOME/data/index
ver=ens92
org=hs
adapter=$idx_dir/oligo/anchored-oligo.fa
bed_exon=$idx_dir/$ver/${org}.region.exon.bed
bed_gene=$idx_dir/$ver/${org}.region.gene.bed
bed_intergenic=$idx_dir/$ver/${org}.region.intergenic.bed
bed_intron=$idx_dir/$ver/${org}.region.intron.bed
genome_size=$idx_dir/$ver/${org}.genome_size
gtf=$idx_dir/$ver/${org}.gtf
ht2_index=$idx_dir/$ver/ht2/${org}

#-- programs --#
bedtools=bedtools
fastqc=fastqc
featureCounts=featureCounts
flexbar=flexbar
hisat2=hisat2
samtools=samtools
stringtie=stringtie
bam2wig=bam2wig.py

#-- pipeline programs --#
sdir=""
concatenate_qc=${sdir}concatenate-qc.py
region_tagger=${sdir}region-tagger.py
extract_isoform_expression=${sdir}extract-isoform-expression.py
join_expression_matrix=${sdir}join-expression-matrix.py
qc_alignment=${sdir}qc-alignment.py
qc_coverage=${sdir}qc-coverage.py
qc_distribution=${sdir}qc-distribution.py
qc_insert_size=${sdir}qc-insert-size.py
wrap_featureCounts=${sdir}wrap-featureCounts.py

#-- prefix --#
raw_pref=raw
cln_pref=cln
qc_dist_pref=qc_distribution
qc_aln_pref=qc_alignment
qc_isize_pref=qc_insert_size
qc_cov_pref=qc_coverage
aln_pref=aln

#-- file name --#
raw1=${raw_pref}_1.fastq.gz
raw2=${raw_pref}_2.fastq.gz
cln1=${cln_pref}_1.fastq
cln2=${cln_pref}_2.fastq
cln1gz=${cln_pref}_1.fastq.gz
cln2gz=${cln_pref}_2.fastq.gz

fqc_raw1=${raw_pref}_1_fastqc.zip
fqc_raw2=${raw_pref}_2_fastqc.zip
fqc_cln1=${cln_pref}_1_fastqc.zip
fqc_cln2=${cln_pref}_2_fastqc.zip

qc_out_dist1=${qc_dist_pref}_1.txt
qc_out_dist2=${qc_dist_pref}_2.txt
qc_out_aln1=${qc_aln_pref}_1.txt
qc_out_aln2=${qc_aln_pref}_2.txt
qc_out_isize=${qc_isize_pref}.txt
qc_out_cov=${qc_cov_pref}.txt

aln=${aln_pref}.bam
aln_hits=${aln_pref}.hits.bam
aln_gene=${aln_pref}.gene.bam
aln_exon=${aln_pref}.exon.bam
aln_intron=${aln_pref}.intron.bam
aln_intergenic=${aln_pref}.intergenic.bam
aln_ercc=${aln_pref}.ercc.bam
aln_unmap=${aln_pref}_unmapped.bam

gene_abund=gene_abund.tab
iso_abund=isoform_abund.tab
trans_gtf=transcript.gtf
tpm_matrix=tpm_matrix.txt
cnt_matrix=count_matrix.txt

