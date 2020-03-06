#-- summary --#
$join_expression_matrix -o $summary_dir/$tpm_matrix $gtf $smpl
$wrap_featureCounts -T $cpu -o $summary_dir/$cnt_matrix $gtf $smpl
$concatenate_qc -o $summary_dir/${qc_aln_pref}.txt $smpl $qc_out_aln
$concatenate_qc -o $summary_dir/${qc_dist_pref}.txt $smpl $qc_out_dist
$concatenate_qc -o $summary_dir/${qc_out_cov} $smpl $qc_out_cov

