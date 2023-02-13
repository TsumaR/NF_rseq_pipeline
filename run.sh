#!/bin/bash
#$ -cwd
#$ -e log
#$ -o log
#$ -pe def_slot 1
#$ -l s_vmem=128G,mem_req=128G

source ~/.bash_profile

export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -Xmx100g"
export TOWER_ACCESS_TOKEN=YOUR_ACCESS_TOKEN
export TOWER_ACCESS_TOKEN=eyJ0aWQiOiA2NzIwfS5mMjI4OGQ5MWEwODE1NDM5ZTk1ZTQzYTdhNTVlOWMwNGJlZmNjMjMy
export NXF_VER=22.04.3

nextflow run nextflow/main.nf -c run.config -profile main -resume -with-report -with-timeline -with-tower
nextflow run nextflow/hisat2.nf -c run.config -profile main -resume -with-report -with-timeline -with-tower
nextflow run nextflow/stringtie.nf -c run.config -profile main -resume -with-report -with-timeline -with-tower
nextflow run nextflow/qc_for_R.nf -c run.config -profile main -resume -with-report -with-timeline -with-tower
nextflow run nextflow/summary.nf -c run.config -profile main -resume -with-report -with-timeline -with-tower
