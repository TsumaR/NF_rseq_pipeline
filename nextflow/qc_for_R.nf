#!/home/myne812/bin/nextflow

proj_id = params.project_id
fastq_filelist = params.fastq_filelist

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

Channel
    .fromPath(fastq_filelist)
    .splitCsv(header: true, sep: "\t")
    .map{ 
        it.Sample_ID
    }
    .into{Sample_id; Sample_ids;}

bam_files = Channel
    .fromPath("output_" + params.project_id + "/**/04_tagging/*_trim.sort.bam")
    .map { [it.baseName.replaceAll('_R1', '').replaceAll('_trim.sort', ''), file(file(it).parent.toString().replaceAll('/04_tagging','')).name, it]}
    // {sample_ids; run_ids; file_path}

fastqc_files = Channel
    .fromFilePairs("./output_" + params.project_id + "/**/*_{R1,R2}*fastqc.zip", size:4, flat : true)

input_paths = bam_files
    .join(fastqc_files)

input_paths
    .into{
        qc_alignment_input
        qc_distribution_input
        qc_insert_input
    }

    
process qc_alignment {

    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set sample_id, run_ids, file(bam), file(cln_1), file(raw_1), file(cln_2), file(raw_2) from qc_alignment_input
    path qc_alignment from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-alignment.py"

    output:
    file("*.txt")

    script:
    """
    python $qc_alignment -o ${sample_id}_R1_alignment.txt $raw_1 $cln_1 $bam
    python $qc_alignment --read2 -o ${sample_id}_R2_alignment.txt $raw_1 $cln_1 $bam
    """
}

process qc_distribution {
    tag{"${run_ids}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set sample_id, run_ids, file(bam), file(cln_1), file(raw_1), file(cln_2), file(raw_2) from qc_distribution_input
    path qc_distribution from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-distribution.py"

    output:
    file("*.txt")

    script:
    """
    python $qc_distribution -o ${sample_id}_R1_distribution.txt $cln_1 $bam
    python $qc_distribution --read2 -o ${sample_id}_R2_distribution.txt $cln_2 $bam
    """
}

process qc_insert_size {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set sample_id, run_ids, file(bam), file(cln_1), file(raw_1), file(cln_2), file(raw_2) from qc_insert_input
    path qc_insert_size from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-insert-size.py"

    output:
    file("*txt")

    script:
    """
    python $qc_insert_size -o ${sample_id}_insertsize.txt $bam
    """
}


isoform_files = Channel
    .fromPath("output_" + params.project_id + "/**/07_ex_isoform/*tab")
    .map { [it.baseName.replaceAll('_R1_isoform_abund', ''), file(file(it).parent.toString().replaceAll('/07_ex_isoform','')).name, it]}

wig_files = Channel
    .fromPath("output_" + params.project_id + "/**/05_bam2wig/*bw")
    .map { [it.baseName.replaceAll('_R1', ''), it]}

gtf_file = Channel
    .from(params.gtf_file)
    .map{ [it[0], file(it[1])] }

qc_coverage_input = isoform_files
    .join(wig_files)
    .combine(gtf_file)

process qc_coverage {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.3"

    input:
    val proj_id
    set sample_id, run_ids, file(iso), file(bw), gtf_name, file(gtf) from qc_coverage_input
    path qc_coverage from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-coverage.py"

    output:
    file("*.txt")

    script:
    """
    python $qc_coverage --isoform $iso -o ${sample_id}_coverage.txt $gtf $bw
    """
}

