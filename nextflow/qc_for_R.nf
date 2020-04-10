#!/home/myne812/bin/ nextflow

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
        //run_id, sample_id
    }
    .into{Sample_id; Sample_ids;}

Channel
    .fromFilePairs("./output_"+params.project_id + "/**/S*_{R1,R2}*fastqc.zip", size:4, flat:true) 
    //.toSortedList( {a, b -> a[1] <=> b[1]} )
    .into{cln_qc_files; cln_qc_files_;}
Channel
    .fromPath("output_" + params.project_id + "/**/04_tagging/", type: "dir")
    .into{bam_dir; bam_dir_;}

qc_alignment_input = bam_dir
    .combine(cln_qc_files)


process qc_alignment {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set bam_dir, sample_id, file(cln_1), file(raw_1), file(cln_2), file(raw_2) from qc_alignment_input
    path qc_alignment from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-alignment.py"

    output:
    set bam_dir, sample_id, file(cln_1), file(raw_1), file(cln_2), file(raw_2) into alignment_output
    file("*.txt")

    script:
    """
    python $qc_alignment -o ${sample_id}_R1_alignment.txt $raw_1 $cln_1 ${bam_dir}/${sample_id}_R1_trim.sort.bam
    python $qc_alignment --read2 -o ${sample_id}_R2_alignment.txt $raw_1 $cln_1 ${bam_dir}/${sample_id}_R1_trim.sort.bam
    """
}

process qc_distribution {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set bam_dir, sample_id, cln_1, raw_1, cln_2, raw_2 from alignment_output
    path qc_distribution from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-distribution.py"

    output:
    set sample_id, bam_dir into distribution_output
    file("*.txt")

    script:
    """
    python $qc_distribution -o ${sample_id}_R1_distribution.txt $cln_1 ${bam_dir}/${sample_id}_R1_trim.sort.bam
    python $qc_distribution --read2 -o ${sample_id}_R2_distribution.txt $cln_2 ${bam_dir}/${sample_id}_R1_trim.sort.bam
    """
}

process qc_insert_size {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set sample_id, bam_dir from distribution_output
    path qc_insert_size from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-insert-size.py"

    output:
    val sample_id into insert_size_output
    file("*txt")

    script:
    """
    python $qc_insert_size -o ${sample_id}_insertsize.txt ${bam_dir}/${sample_id}_R1_trim.sort.bam
    """
}

Channel
    .fromPath("output_" + params.project_id + "/**/07_ex_isoform/", type: "dir")
    .set{iso_dir}

Channel
    .fromPath("output_" + params.project_id + "/**/05_bam2wig/", type: "dir")
    .set{bw_dir}

gtf_file = Channel
    .from(params.gtf_file)
    .map{ [it[0], file(it[1])] }

iso_abundant_input = gtf_file
    .combine(gtf_file)
    .combine(iso_dir)
    .combine(bw_dir)


process qc_coverage {
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/08_qc", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.3"

    input:
    val proj_id
    set run_ids, sample_id, gtf_name, file(gtf), iso_dir, bw_dir from iso_abundant_input
    path qc_coverage from workflow.scriptFile.parent.parent + "/rnaseq-pipeline-4.0.2/qc-coverage.py"

    output:
    file("*.txt")

    script:
    """
    python $qc_coverage --isoform ${iso_dir}/${sample_id}*.tab -o ${sample_id}_coverage.txt $gtf ${bw_dir}/${sample_id}*bw
    """
}