#!/home/myne812/bin/ nextflow

proj_id = params.project_id
fastq_filelist = params.fastq_filelist

Channel
    .from(params.run_ids)
    .map{[it[0]]}
    .into{run_ids; run_ids_;}

Channel
    .from(params.gtf_file)
    .map{ [it[0], file(it[1])] }
    .into{gtf_file; gtf_file_;}

Channel
    .fromPath(params.fastq_filelist)
    .into{sample_list; sample_list_;}

Channel
    .fromPath("output_" + params.project_id + "/**/06_stringtie/", type: "dir")
    .into{tab_dir; tab_dir_;}

gtf_file
    .combine(sample_list)
    .combine(tab_dir)
    .into{joint_conditions; joint_conditions_;}

process joint_ex_mx {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set gtf_name, file(gtf), file(smpl), tab_dir from joint_conditions
    path joint_expression_matrix from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.3-nextflow/joint-expression-matrix.py"

    output:
    file("tpm_matrix.txt")

    script:
    """
    python $joint_expression_matrix -o tpm_matrix.txt $gtf $smpl -i $tab_dir
    """
}

Channel
    .fromPath("output_" + params.project_id + "/**/04_tagging/", type: "dir")
    .into{bam_dir; bam_dir_;}
cpu = Channel
    .from(params.cpu)
    .map{ [it[0]] }

gtf_file_
    .combine(sample_list_)
    .combine(cpu)
    .combine(bam_dir)
    .into{wrap_conditions; wrap_conditions_;}

process wrap_count {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/subread:1.0"

    input:
    val proj_id
    set gtf_name, file(gtf), file(smpl), cpu, bam_dir from wrap_conditions
    path wrap_featureCounts from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.3-nextflow/wrap-featureCounts.py"

    output:
    set file(smpl), cpu into wrap_count_output
    file("count_matrix.txt")

    script:
    """
    python $wrap_featureCounts -p -T $cpu -o count_matrix.txt $gtf $smpl -i $bam_dir
    """
}

Channel
    .fromPath("output_" + params.project_id + "/**/08_qc/", type: "dir")
    .into{qc_out_dir; qc_out_dir_;}

wrap_count_output
    .combine(qc_out_dir)
    .into{concentrate_qc_input; bef;}

process concatenate_qc {
    tag{"${proj_id}"}
    publishDir "./summary", mode: 'copy', overwrite: true

    container "tsumar/rseq_python:1.2"

    input:
    val proj_id
    set file(smpl), cpu, qc_dir from concentrate_qc_input
    path concatenate_qc from workflow.scriptFile.parent.parent + "/rseq-pipe-4.0.3-nextflow/concatenate-qc.py"

    output:
    file("*.txt")

    script:
    """
    python $concatenate_qc -t alignment -o qc_alignment.txt -pe LR $smpl $qc_dir
    python $concatenate_qc -t distribution -o qc_distribution.txt -pe LR $smpl $qc_dir
    python $concatenate_qc -t insertsize -o qc_insertsize.txt $smpl $qc_dir
    python $concatenate_qc -t coverage -o qc_coverage.txt $smpl $qc_dir
    """
}

