#!/bin/bash
#$ -cwd
#$ -e log
#$ -o log
#$ -pe def_slot 4
#$ -l s_vmem=8G,mem_req=8G

module load java/8
source ~/.bash_profile

PWD=`pwd`
fq_dir=${PWD}/fastq
smpl_lst=(`cat ${PWD}/ifiles.txt|xargs`)
# echo ${smpl_lst[@]}

echo -e "Sample_ID\tSequence_ID\tFastq1\tFastq2" > sample.txt
for s_name in ${smpl_lst[@]}; do
  echo $s_name
  # fqs=(`find ${fq_dir} -name "${s_name}*"`)
  r1=(`ls ${fq_dir}/${s_name}*R1*`)
  r2=(`ls ${fq_dir}/${s_name}*R2*`)
  echo $r1
  echo $r2

  # R1=${fqs[0]}
  # R2=${fqs[1]}
  echo -e "$s_name\t"NF_rseq_pipeline"\t$r1\t$r2" >> sample.txt
done
