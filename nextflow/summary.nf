#!/home/myne812/bin/ nextflow

proj_id = params.project_id
fastq_filelist = params.fastq_filelist

// Channel
    // .from(params.run_ids)
    // .map{[it[0]]}
    // .into{run_ids; run_ids_;}

Channel
    .from(params.gtf_file)
    .map{ [it[0], file(it[1])] }
    .into{gtf_file; gtf_file_;}

Channel
    .fromPath(params.fastq_filelist)
    .into{sample_list; sample_list_; sample_list_ad}

Channel
    .fromPath("output_" + params.project_id, type:"dir")
    .into{out_dir_; out_dir_second; out_dir_third;}

joint_conditions = gtf_file
    .combine(sample_list)
    .combine(out_dir_)

process joint_ex_mx {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set gtf_name, file(gtf), file(smpl), out_dir from joint_conditions
    path joint_expression_matrix from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.2-nextflow/joint-expression-matrix.py"

    output:
    file("tpm_matrix.txt")

    script:
    """
    python $joint_expression_matrix -o tpm_matrix.txt $gtf $smpl -i $out_dir
    """
}

cpu = Channel
    .from(params.cpu)
    .map{ [it[0]] }

wrap_count_input = gtf_file_
    .combine(sample_list_)
    .combine(cpu)
    .combine(out_dir_second)

process wrap_count {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/subread:1.0"

    input:
    val proj_id
    set gtf_name, file(gtf), file(smpl), cpu, out_dir from wrap_count_input
    path wrap_featureCounts from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.2-nextflow/wrap-featureCounts.py"

    output:
    file("count_matrix.txt")

    script:
    """
    python $wrap_featureCounts -p -T $cpu -o count_matrix.txt $gtf $smpl -i $out_dir
    """
}

concatenate_input = sample_list_ad
    .combine(out_dir_third)

process concatenate_qc {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set file(smpl), out_dir from concatenate_input
    path concatenate_qc from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.2-nextflow/concatenate-qc.py"

    output:
    file("*.txt")

    script:
    """
    python $concatenate_qc -t alignment -o qc_alignment.txt -pe LR $smpl $out_dir
    python $concatenate_qc -t distribution -o qc_distribution.txt -pe LR $smpl $out_dir
    python $concatenate_qc -t insertsize -o qc_insertsize.txt $smpl $out_dir
    python $concatenate_qc -t coverage -o qc_coverage.txt $smpl $out_dir
    """
}

