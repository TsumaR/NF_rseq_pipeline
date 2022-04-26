#!/usr/bin/env nextflow

proj_id = params.project_id
fastq_filelist = params.fastq_filelist

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

bam_files = Channel
    .fromPath("output_" + params.project_id + "/**/04_tagging/*_trim.sort.bam")
    .map { [file(file(it).parent.toString().replaceAll('/04_tagging','')).name, it.baseName.replaceAll('_trim', '').replaceAll('.sort', ''), it]}


bam_files
    .into{
        bam_files_test
        bam_files_input
        bam_files_to_count
    }

gtf_file = Channel
    .from(params.gtf_file)
    .map{ [it[0], file(it[1])] }

cpu = Channel
    .from(params.cpu)
    .map{ [it[0]] }

stringtie_conditions = bam_files_input
    .combine(gtf_file)
    .combine(cpu)

process run_stringtie {
 
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/06_stringtie", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/stringtie:2.1.7--h978d192_0"

    input:
    val proj_id
    set run_ids, bam_name, file(bam), gtf_name, file(gtf), cpu from stringtie_conditions

    output:
    set run_ids, bam_name, file(bam), file("${bam_name}_abund.tab"), file("${bam_name}.gtf") into stringtie_output
    file("${bam_name}_abund.tab")
    file("${bam_name}.gtf")

    
    script:
    """
    stringtie -e -p $cpu -A ${bam_name}_abund.tab -o ${bam_name}.gtf -G $gtf $bam
    """
}

process extract_isoform {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/07_ex_isoform", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.1"
    input:
    val proj_id
    set run_ids, bam_name, file(bam), abund_tab, file(trans_gtf) from stringtie_output
    path extract_isoform_expression from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/extract-isoform-expression.py"

    output:
    set run_ids, file("${bam_name}_isoform_abund.tab") into extract_output
    file("${bam_name}_isoform_abund.tab")

    script:
    """
    python $extract_isoform_expression $trans_gtf > ${bam_name}_isoform_abund.tab
    """
}

/*
bam2wig_ch = extract_output
    .combine(hisat2_index)

process bam2wig {

    tag {"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/02_hisat2", mode: 'copy', overwrite: true

    input:
    val project_id
    set run_ids, bam_name, bam_file, abund_tab, genome_size  

    script:
    """
    bam2wig.py -i $bam_file -s $genome_size -o $bam_name
    """
}
*/
