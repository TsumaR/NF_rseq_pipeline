#!/bin/bash
#$ -cwd
#$ -e log
#$ -o log
#$ -pe def_slot 4
#$ -l s_vmem=8G,mem_req=8G

~/bin/nextflow run nextflow/main.nf -c run.config -resume -with-report log.01.main.html
~/bin/nextflow run nextflow/hisat2.nf -c run.config -resume -with-report log.02.hisat2.html
~/bin/nextflow run nextflow/stringtie.nf -c run.config -resume -with-report log.03.stringtie.html
~/bin/nextflow run nextflow/qc_for_R.nf -c run.config -resume -with-report log.04.qc_for_R.html
~/bin/nextflow run nextflow/summary.nf -c run.config -resume -with-report log.05.summary.html
