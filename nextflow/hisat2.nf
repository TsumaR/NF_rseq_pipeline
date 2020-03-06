#!/home/myne812/bin/ nextflow

proj_id = params.project_id
fastq_filelist = params.fastq_filelist

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

fastq_files = Channel
    .fromFilePairs("./output_" + params.project_id + "/**/*_trim_{1,2}.fastq.gz", flat: true) { file -> file.name.replaceAll(/_1|_2/,'').replaceAll('_trim', '').replaceAll('.fastq.gz', '') }
    .map{[file(file(it[1]).parent.toString().replaceAll('/01_flexbar','')).name, it[0], it[1], it[2]]}

fastq_files
    .into{
        fastq_files_test
        fastq_files_input
        fastq_files_to_count
    }

hisat2_index = Channel
    .from(params.hisat2_index)
    .map{ [it[0], it[1], file(it[1]+"*")] }
cpu = Channel
    .from(params.cpu)
    .map{ [it[0]] }
hisat2_conditions = fastq_files_input
    .combine(hisat2_index)
    .combine(cpu)

/*
withName ：foo  { 
        container  =  'image_name_1' 
    } 
    withName ：bar  { 
        container  =  'image_name_2' 
    } 
} 
*/

process run_hisat2 {

    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/03_hisat2", mode: 'copy', overwrite: true
    //container "quay.io/biocontainers/hisat2:2.1.0--py37hc9558a2_4"
    container "kathrinklee/rna-seq-pipeline-hisat2:latest"

    input:
    val proj_id
    set run_ids, fastq_name, file(fastq_L), file(fastq_R), index_name, path(index), file(index_files), cpu_number from hisat2_conditions

    output:
    set run_ids, fastq_name, file("*.bam"), cpu_number into hisat2_output
    file "*.bam" 

    script:
    """
    hisat2 -p $cpu_number -x $index -1 $fastq_L -2 $fastq_R | samtools view -b | samtools sort -@ $cpu_number -o ${fastq_name}_trim.sort.bam
    """ 
}

gene_file = Channel
    .from(params.gene_file)
    .map{ [it[0], file(it[1])] }

exon_file = Channel
    .from(params.exon_file)
    .map{ [it[0], file(it[1])] }

tag_opt = hisat2_output
    .combine(gene_file)
    .combine(exon_file)

process tagging {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/04_tagging", mode: 'copy', overwrite: true
    container "tsumar/sambedtools:1.1"

    input:
    val proj_id
    set run_ids,fastq_name, file(aln), cpu, bed_gene, file(gene_file), bed_exon, file(exon_file) from tag_opt
    path region_tagger from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/region-tagger.py"

    output:
    set run_ids, fastq_name, file("${aln}"), file("${aln}.bai") into tagging_output
    file("${aln}")

    script:
    """
    samtools view -b -f4 $aln > ${fastq_name}_unmapped.bam
    samtools view -h -F4 $aln | grep -v '\t'ERCC | samtools view -b > ${fastq_name}.hits.bam
    samtools view -h -F4 $aln | grep -e '^@' -e '\t'ERCC | samtools view -b > ${fastq_name}.ercc.bam
    bedtools intersect -wa -u -a ${fastq_name}.hits.bam -b $gene_file > ${fastq_name}.gene.bam
    bedtools intersect -wa -v -a ${fastq_name}.hits.bam -b $gene_file > ${fastq_name}.intergenic.bam
    bedtools intersect -wa -u -a ${fastq_name}.gene.bam -b $exon_file > ${fastq_name}.exon.bam
    bedtools intersect -wa -v -a ${fastq_name}.gene.bam -b $exon_file > ${fastq_name}.intron.bam
    /opt/conda/bin/python2.7 $region_tagger -o ${fastq_name}.exon.bam ${fastq_name}.exon.bam exon
    /opt/conda/bin/python2.7 $region_tagger -o ${fastq_name}.intron.bam ${fastq_name}.intron.bam intron
    /opt/conda/bin/python2.7 $region_tagger -o ${fastq_name}.intergenic.bam ${fastq_name}.intergenic.bam intergenic
    /opt/conda/bin/python2.7 $region_tagger -o ${fastq_name}.ercc.bam ${fastq_name}.ercc.bam ercc
    samtools merge -f -@ $cpu $aln ${fastq_name}.exon.bam ${fastq_name}.intron.bam ${fastq_name}.intergenic.bam ${fastq_name}.ercc.bam
    samtools sort -@ $cpu -o $aln $aln
    samtools index $aln
    """
}

chr_size = Channel
    .from(params.ref_chrsize)
    .map{ [it[0], file(it[1])] }

bam2wig_ch = tagging_output
    .combine(chr_size)

process bam2wig {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/05_bam2wig", mode: 'copy', overwrite: true

    input:
    val proj_id
    set run_ids, bam_name, file(bam), file(bai_file),chr_size, file(genome_size) from bam2wig_ch
    container "tsumar/bamtobigwig:3.0.0"

    output:
    file("${bam_name}.bw")

    script:
    """
    bam2wig.py -i $bam -s $genome_size -o $bam_name
    """
}