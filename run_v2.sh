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

nextflow run main_v2.nf -profile main -resume -with-report -with-timeline -with-tower -with-dag flowdiagram.mmd
