#!/bin/bash
#$ -cwd
#$ -e log
#$ -o log
#$ -pe def_slot 4
#$ -l s_vmem=8G,mem_req=8G

module load java/8
source ~/.bash_profile

PWD=`pwd`
fq_dir=./fastq
smpl_lst=(`cat ${PWD}/ifiles.txt|xargs`)

echo -e "Sample_ID\tSequence_ID\tFastq1\tFastq2" > sample.txt
for s_name in ${smpl_lst[@]}; do
  fqs=(`find ${fq_dir} -name "${s_name}*"`)
  R1=${fqs[0]}
  R2=${fqs[1]}
  echo -e "$s_name\t"NF_rseq_pipeline"\t$R1\t$R2" >> sample.txt
done
