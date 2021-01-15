# NF_rseq_pipeline 
Nextflowで記載したRNA-seq pipeline。解析に必要なカウントデータなどのテーブルデータを複数作成する。
全てdocker環境で実施するようにしてあるので，Nextflowとdocker(singularity)さえあれば他の環境構築は必要ない。 
東大スパコンとOILのスパコンにはsingularityが標準で存在するので実質必要なのはNextflowだけである。

## Requirements
- [Nextflow](https://www.nextflow.io/)

## Getting started

### 1. Installing Nextflow

1. Make sure 8 or later is installed on your computer by using the command: `java -version`
2. Enter the below commands in your terminal (The command creates a file nextflow in `~/bin`)
この作業は計算機ごとに一度行えば今後行う必要はありません。

```
module load java/8
mkdir -p ~/bin
cd ~/bin
wget -qO- https://get.nextflow.io | bash
``` 
3. Run the classic Hello world by entering the following command: `~/bin/nextflow run hello` 

### 2. singularity

You have to open path of singularity
If you use shirokane, you add below sentense to `.bash_profile`
この作業は計算機ごとに一度行えば今後行う必要はありません。

```
export PATH=/usr/local/package/singularity/3.2.1/bin:$PATH
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity
```

### 3. Clone this repository

`WORKING_DIR`には自分の解析したい実験名のディレクトリを指定してください。

```
cd  "WORKING_DIR"
git clone https://github.com/TsumaR/NF_rseq_pipeline.git
```

### 4. Modifying config file 

You can edit `local.config` file to change the name of the analysis result.
これにより、project nameなどを変更することができますが、必須の作業ではありません。

```
vim local.config
```

### 5. Preparing sample file

### 6. Storing FASTQ file

### 5. Run the pipeline

```
~/bin/nextflow run nextflow/main.nf -c run.config -resume -with-report log.01.main.html
~/bin/nextflow run nextflow/hisat2.nf -c run.config -resume -with-report log.02.hisat2.html
~/bin/nextflow run nextflow/stringtie.nf -c run.config -resume -with-report log.03.stringtie.html
~/bin/nextflow run nextflow/qc_for_R.nf -c run.config -resume -with-report log.04.qc_for_R.html
~/bin/nextflow run nextflow/summary.nf -c run.config -resume -with-report log.05.summary.html
``` 

or

```
qsub run..sh
```

## Information 
Version of packages

2020/04/09，東大スパコンのデフォルトのJavaだと動かないエラー発生

```
#現在のバージョンを確認
$ which java
/usr/local/package/java/10_2018-03-20/bin/java

#javaのバージョン8を読み込み
$ module load java/8
$ which java
/usr/local/package/java/1.8.0_181/bin/java
``` 
