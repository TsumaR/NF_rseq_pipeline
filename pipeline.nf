#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process run_flexbar {

  publishDir path: 'result/flexbar', mode: 'copy', overwrite: true

  input:
    tuple val(sampleid), path(reads), path(adapter)
  output:
    tuple val(sampleid), path("${sampleid}_qc_1.fastq.gz"), path("${sampleid}_qc_2.fastq.gz"), path("${sampleid}_raw_1.fastq.gz"), path("${sampleid}_raw_2.fastq.gz"), emit: tohisat2
    tuple val(sampleid), path("${sampleid}_qc_1.fastq.gz"), path("${sampleid}_qc_2.fastq.gz"), path("${sampleid}_raw_1.fastq.gz"), path("${sampleid}_raw_2.fastq.gz"), emit: tofastqc

  shell:
  """
    echo ${sampleid}
    flexbar \
      -n !{task.cpus} \
      -z GZ -u 10 -at ANY --htrim-right AT --htrim-min-length 10 --htrim-error-rate 0.1 \
      -a ${adapter} \
      -t ${sampleid}_qc \
      -r !{reads[0]} -p !{reads[1]}
    cp !{reads[0]} !{sampleid}_raw_1.fastq.gz
    cp !{reads[1]} !{sampleid}_raw_2.fastq.gz
  """
}

process run_fastqc {

  publishDir path: 'result/fastqc', mode: 'copy', overwrite: true

  input:
    tuple val(sampleid), path(qc1), path(qc2), path(raw1), path(raw2)
  output:
    tuple val(sampleid), path(qc1), path(qc2), path(raw1), path(raw2), path("${sampleid}_qc_1_fastqc.zip"), path("${sampleid}_qc_2_fastqc.zip"), path("${sampleid}_raw_1_fastqc.zip"), path("${sampleid}_raw_2_fastqc.zip"), emit: out_fastqc

  shell:
  '''
    zcat !{qc1} | fastqc stdin:!{sampleid}_qc_1
    zcat !{qc2} | fastqc stdin:!{sampleid}_qc_2
    zcat !{raw1} | fastqc stdin:!{sampleid}_raw_1
    zcat !{raw2} | fastqc stdin:!{sampleid}_raw_2
  '''
}

process run_hisat2 {

  publishDir path: 'result/hisat2', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(qc1), path(qc2), path(raw1), path(raw2), path(qc1_repo), path(qc2_repo), path(raw1_repo), path(raw2_repo)
  output:
    tuple val(sampleid), path("${sampleid}.bam"), path(qc1), path(qc2), path(raw1), path(raw2), path(qc1_repo), path(qc2_repo), path(raw1_repo), path(raw2_repo), emit: totagging

  shell:
  '''
    hisat2 -p !{task.cpus} -x !{params.hisat2_index} \
      -1 !{qc1} -2 !{qc2} | samtools view -b | samtools sort -@ !{task.cpus} > !{sampleid}.bam
  '''
}

process tagging {

  publishDir path: 'result/tagging', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(hisat2bam), path(qc1), path(qc2), path(raw1), path(raw2), path(qc1_repo), path(qc2_repo), path(raw1_repo), path(raw2_repo), path(gene_file), path(exon_file), path(regiontagger)

  output:
    tuple val(sampleid), path("${hisat2bam}"), path("${sampleid}.bam.bai"), path(qc1), path(qc2), path(raw1), path(raw2),emit: tobam2wig
    tuple val(sampleid), path("${hisat2bam}"), path("${sampleid}.bam.bai"), path(qc1_repo), path(qc2_repo), path(raw1_repo), path(raw2_repo), emit: toqc_alignment

  shell:
  '''
    samtools view -b -f4 !{hisat2bam} > !{sampleid}.unmapped.bam
    samtools view -h -F4 !{hisat2bam} | grep -v '\t'ERCC | samtools view -b > !{sampleid}.hits.bam
    samtools view -h -F4 !{hisat2bam} | grep -e '^@' -e '\t'ERCC | samtools view -b > !{sampleid}.ercc.bam
    bedtools intersect -wa -u -a !{sampleid}.hits.bam -b !{gene_file} > !{sampleid}.gene.bam
    bedtools intersect -wa -v -a !{sampleid}.hits.bam -b !{gene_file} > !{sampleid}.intergenic.bam
    bedtools intersect -wa -u -a !{sampleid}.gene.bam -b !{exon_file} > !{sampleid}.exon.bam
    bedtools intersect -wa -v -a !{sampleid}.gene.bam -b !{exon_file} > !{sampleid}.intron.bam
    /opt/conda/bin/python2.7 !{regiontagger} -o !{sampleid}.exon.bam !{sampleid}.exon.bam exon
    /opt/conda/bin/python2.7 !{regiontagger} -o !{sampleid}.intron.bam !{sampleid}.intron.bam intron
    /opt/conda/bin/python2.7 !{regiontagger} -o !{sampleid}.intergenic.bam !{sampleid}.intergenic.bam intergenic
    /opt/conda/bin/python2.7 !{regiontagger} -o !{sampleid}.ercc.bam !{sampleid}.ercc.bam ercc
    samtools merge -f -@ !{task.cpus} !{hisat2bam} !{sampleid}.exon.bam !{sampleid}.intron.bam !{sampleid}.intergenic.bam !{sampleid}.ercc.bam
    samtools sort -@ !{task.cpus} -o !{hisat2bam} !{hisat2bam}
    samtools index !{hisat2bam}
  '''
}

process bam2wig {
  publishDir path: 'result/bam2wig', mode: 'copy', overwrite: true
  input:
    tuple  val(sampleid), path(alnbam), path(bambai), path(qc1), path(qc2), path(raw1), path(raw2), path(genome_size)
  output:
    tuple val(sampleid), path("${sampleid}.bam"), path("${sampleid}.bam.bai"), path("${sampleid}.bw"), path(qc1), path(qc2), path(raw1), path(raw2),emit: tostringtie
  shell:
  '''
    bam2wig.py -i !{alnbam} -s !{genome_size} -o !{sampleid}
  '''

}

process run_stringtie {
  publishDir path: 'result/stringtie', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(alnbam), path(bambai), path(bw),path(qc1), path(qc2), path(raw1), path(raw2),path(gtf)
  output:
    tuple val(sampleid), path("${sampleid}.gtf"), path("${sampleid}_abund.tab"), path(bw), emit: toisoform
    //tuple val(sampleid), path(alnbam), path(qc1), path(qc2), path(raw1), path(raw2), emit: toalignment

  shell:
  '''
    stringtie -e -p !{task.cpus} -A !{sampleid}_abund.tab -o !{sampleid}.gtf -G !{gtf} !{alnbam}
  '''
}

process run_extract_isoform {
  publishDir path: 'result/isoform', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(smplgtf), path(smpl_abund), path(bw), path(Pyextract_isoform)
  output:
    tuple val(sampleid), path("${sampleid}_isoform_abund.tab"), path(bw), emit: toqc_coverage

  shell:
  '''
    python !{Pyextract_isoform} !{smplgtf} > !{sampleid}_isoform_abund.tab   
  '''
}


process qc_alignment {
  publishDir path: 'result/qc_forR', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(bam), path(bai), path(cln1), path(cln2), path(raw1), path(raw2), path(Pyqc_alignment)
  output:
    tuple val(sampleid), path("${sampleid}_R1_alignment.txt"), path("${sampleid}_R2_alignment.txt")
    tuple val(sampleid), path(cln1), path(cln2), path(raw1), path(raw2), path(bam), emit: todistribution

  shell:
  '''
    python !{Pyqc_alignment} -o !{sampleid}_R1_alignment.txt !{raw1} !{cln1} !{bam}
    python !{Pyqc_alignment} --read2 -o !{sampleid}_R2_alignment.txt !{raw1} !{cln1} !{bam}
  '''
}

process qc_distribution {
  publishDir path: 'result/qc_forR', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(cln1), path(cln2), path(raw1), path(raw2), path(bam), path(Pyqc_distribution)
  output:
    tuple val(sampleid), path("${sampleid}_R1_distribution.txt"), path("${sampleid}_R2_distribution.txt")
    tuple val(sampleid), path(bam), emit: toinsert_size
  shell:
  '''
    python !{Pyqc_distribution} -o !{sampleid}_R1_distribution.txt !{cln1} !{bam}
    python !{Pyqc_distribution} --read2 -o !{sampleid}_R2_distribution.txt !{cln2} !{bam}
  '''
}

process qc_insert_size {
  publishDir path: 'result/qc_forR', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(bam), path(Pyqc_insert_size)
  output:
    tuple val(sampleid), path("${sampleid}_insertsize.txt")
  shell:
  '''
    python !{Pyqc_insert_size} -o !{sampleid}_insertsize.txt !{bam}
  '''
}

process qc_coverage {
  publishDir path: 'result/qc_forR', mode: 'copy', overwrite: true
  input:
    tuple val(sampleid), path(iso), path(bw), path(gtf), path(Pyqc_coverage)
  output:
    tuple val(sampleid), path("${sampleid}_coverage.txt")
    tuple val(sampleid), emit: smplid_lists
  shell:
  '''
    python !{Pyqc_coverage} --isoform !{iso} -o !{sampleid}_coverage.txt !{gtf} !{bw}
  '''
}

process joint_ex_mx {
  publishDir path: 'summary', mode: 'copy', overwrite: true
  input:
    tuple path(Pyjoint_ex_mx), path(smpl), path(gtf), path(outdir)
    val(sampleids)

  output:
    tuple path("tpm_matrix.txt")
    tuple path(smpl), path(gtf), path(outdir), emit: smpl_gtf
  shell:
  '''
    python !{Pyjoint_ex_mx} -o tpm_matrix.txt !{gtf} !{smpl} -i !{outdir}
  '''
}

process wrap_count {
  publishDir path: 'summary', mode: 'copy', overwrite: true
  input:
    tuple path(Pywrap_featureCounts), path(smpl), path(gtf), path(outdir)
  output:
    tuple path("count_matrix.txt")
    tuple path(smpl), path(gtf), path(outdir), emit: smpl_gtf
  shell:
  '''
    python !{Pywrap_featureCounts} -p -T !{task.cpus} -o count_matrix.txt !{gtf} !{smpl} -i !{outdir}
  '''
}

process concatnate_qc {
  publishDir path: 'summary', mode: 'copy', overwrite: true
  input:
    tuple path(Pyconcat_qc), path(smpl), path(gtf), path(outdir)
  output:
    tuple path("qc_alignment.txt"), path("qc_distribution.txt"), path("qc_insertsize.txt"), path("qc_coverage.txt")
  shell:
  '''
    python !{Pyconcat_qc} -t alignment -o qc_alignment.txt -pe LR !{smpl} !{outdir}
    python !{Pyconcat_qc} -t distribution -o qc_distribution.txt -pe LR !{smpl} !{outdir}
    python !{Pyconcat_qc} -t insertsize -o qc_insertsize.txt !{smpl} !{outdir}
    python !{Pyconcat_qc} -t coverage -o qc_coverage.txt !{smpl} !{outdir}
  '''
}

workflow {
  //step1
  flexbar_channel = Channel
    .fromFilePairs('fastq/*_R{1,2}_001.fastq.gz')
    .combine(Channel.fromPath(params.adapter))

  run_flexbar(flexbar_channel)
  
  //step2
  run_fastqc(run_flexbar.out.tofastqc)

  //step3
  run_hisat2(run_fastqc.out.out_fastqc)

  //step4
  tagging_input = run_hisat2.out.totagging
    .combine(Channel.fromPath(params.gene_file))
    .combine(Channel.fromPath(params.exon_file))
    .combine(Channel.fromPath(params.regiontagger))
  tagging(tagging_input)

  //step5
  bam2wig_input = tagging.out.tobam2wig
    .combine(Channel.fromPath(params.genome_size))
  bam2wig(bam2wig_input)

  //step6
  stringtie_input = bam2wig.out.tostringtie
    .combine(Channel.fromPath(params.gtf_file))
  run_stringtie(stringtie_input)

  //step7
  isoform_input = run_stringtie.out.toisoform
    .combine(Channel.fromPath(params.extract_isoform))
  run_extract_isoform(isoform_input)

  //step8
  alignment_input = tagging.out.toqc_alignment
    .combine(Channel.fromPath(params.qc_alignment))
  qc_alignment(alignment_input)

  //step9
  distribution_input = qc_alignment.out.todistribution
    .combine(Channel.fromPath(params.qc_distribution))
  qc_distribution(distribution_input)

  //step10
  insert_size_input = qc_distribution.out.toinsert_size
    .combine(Channel.fromPath(params.qc_insert_size))
  qc_insert_size(insert_size_input)

  //step11
  insert_coverage_input = run_extract_isoform.out.toqc_coverage
    .combine(Channel.fromPath(params.gtf_file))
    .combine(Channel.fromPath(params.qc_coverage))
  qc_coverage(insert_coverage_input)

  //step12
  joint_ex_input = Channel
    .fromPath(params.joint_ex_mx)
    .combine(Channel.fromPath('sample.txt'))
    .combine(Channel.fromPath(params.gtf_file))
    .combine(Channel.fromPath('result'))
  joint_ex_mx(joint_ex_input, qc_coverage.out.smplid_lists.collect())

  //step13
  wrap_count_input = Channel
    .fromPath(params.wrap_count)
    .combine(joint_ex_mx.out.smpl_gtf)
  wrap_count(wrap_count_input)

  //step14
  qc_concat_input = Channel
    .fromPath(params.concatnate_qc)
    .combine(wrap_count.out.smpl_gtf)
  concatnate_qc(qc_concat_input)

}
