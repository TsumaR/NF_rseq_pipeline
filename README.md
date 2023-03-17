# NF_rseq_pipeline 
Nextflowで記載したRNA-seq pipeline。解析に必要なカウントデータなどのテーブルデータを複数作成する。
全てdocker環境で実施するようにしてあるので，Nextflowとdocker(singularity)さえあれば他の環境構築は必要ない。 
東大スパコンとOILのスパコンにはsingularityが標準で存在するので実質必要なのはNextflowだけである。

## Requirements
- [Nextflow](https://www.nextflow.io/)
- [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

## run
qsub run.sh
