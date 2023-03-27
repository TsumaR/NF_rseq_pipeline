#!/bin/bash
#$ -cwd
#$ -e log
#$ -o log
#$ -pe def_slot 1
#$ -l s_vmem=128G,mem_req=128G
#$ -M s6z5s8s8h9i5o7c4@htl-org.slack.com

source ~/.bash_profile

export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -Xmx100g"
export TOWER_ACCESS_TOKEN=YOUR_ACCESS_TOKEN
export NXF_VER=22.04.3

# nextflow run main_v2.nf -profile main -resume -with-report -with-timeline -with-tower -with-dag flowdiagram.mmd

# create smpl.txt
echo "Sample_ID" > ./sample.txt
ls fastq/*R1* | sed -e "s/_R1.*//g" | sed -e "s/fastq\///g" >> ./sample.txt

# run
nextflow run pipeline.nf -profile pipeline -resume -with-report -with-timeline -with-tower -with-dag flowdiagram.mmd
