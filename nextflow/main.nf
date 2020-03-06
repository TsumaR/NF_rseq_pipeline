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
        [file(it.Fastq1).parent.toString().split('/')[-2], it.Sample_ID + "_R1", it.Sample_ID + "_R2", file(it.Fastq1), file(it.Fastq2), file(it.Fastq1).baseName.replaceAll('.fastq',''), file(it.Fastq2).baseName.replaceAll('.fastq','')]
        //run_id, fastq_L_name, fastq_R_name, fastq_L, fastq_R, fastq_L_basename, fastq_R_basename
    }
    .into{fastq_files; fastq_files_tmp;}

flexbar_options = Channel
    .from(params.flexbar_condition)
    .map{ [it[0]] }
adapters = Channel
    .from(params.adapter)
    .map{ [it[0], file(it[1])] }
//cpu = Channel
//    .from(params.flexbar_cpu)
//    .map{ [it[0]] }    

fastq_files
    .combine(flexbar_options)
    .combine(adapters)
    //.combine(cpu)
    .into{flexbar_bef_conditions; flexbar_bef_conditions_;}

process run_flexbar {
  
    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/01_flexbar", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/flexbar:3.5.0--hb7ba0dd_3"
  
    input:
    val proj_id
    set run_ids, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), fastq_L_basename, fastq_R_basename, option_name, adapter_name, file(adapterfile) from flexbar_bef_conditions 

    output:
    set run_ids, fastq_L_name, fastq_R_name, file(fastq_L), file(fastq_R), file("${fastq_L_name}_trim_1.fastq.gz"), file("${fastq_L_name}_trim_2.fastq.gz") into flexbar_output

    script:
    """
    flexbar $option_name -a $adapterfile -t ${fastq_L_name}_trim -r $fastq_L -p $fastq_R
    """
}
//flexbar_output.println()


process run_fastqc {

    tag{"${proj_id}"}
    publishDir "output_${proj_id}/${run_ids}/02_fastqc", mode: 'copy', overwrite: true

    container "quay.io/biocontainers/fastqc:0.11.9--0"

    input:
    val proj_id
    set run_ids, fqL_name, fqR_name, file(fastq_L), file(fastq_R), file(cln_L), file(cln_R) from flexbar_output

    output:
    file("*fastqc*")

    script:
    """
    zcat $fastq_L | fastqc stdin:${fqL_name}
    zcat $fastq_R | fastqc stdin:${fqR_name}
    zcat $cln_L | fastqc stdin:${fqL_name}_cln
    zcat $cln_R | fastqc stdin:${fqR_name}_cln
    """
}
