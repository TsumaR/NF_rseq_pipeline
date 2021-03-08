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

デフォルトではマウスに対するパイプラインになっています。
もしヒトサンプルに対してパイプラインを実行したい場合は`run.config`ファイルを変更する必要があります。

```
vim run.config
```

で`run.config`ファイルの編集画面に入り、

```
process.executor = 'uge'
process.memory = '24G'

includeConfig "./config_file/mm.config"
includeConfig "./config_file/pipeline.config"
```
を
```
process.executor = 'uge'
process.memory = '24G'

includeConfig "./config_file/hs.config"
includeConfig "./config_file/pipeline.config"
```
に変更してください。
mmの部分をhsに変更するだけです。

### 5. Preparing Fastq files

最初に、scpコマンドを利用して自分のPC上のFastqファイルを計算機に送ってください
この際、fastqファイルは計算機上で`fastq`ディレクトリに保存するようにしてください。
fastqファイルは`sample名[自由な文字列]_R1_[自由な文字列].fastq.gz`のようにしてください。

```
mkdir fastq
scp [自分のPC上のfastqファイル] ユーザー名@IPアドレス:[実行ディレクトリ]/fastq
```

### 6. Making sample.txt file

サンプルファイルを作成します。
まず、解析に利用したいサンプル名を縦に並べた下のような`ifiles.txt`というファイルを作成してください。
この際、sample名はfastqファイルの先頭部分と等しくなるように設定してください。

このファイルはローカルPCで作成してから`scp`コマンドで送付してもいいですし、`vim`コマンドでスパコン上に作成しても構いません。

```
sample1
sample2
sample3
sample4
```

このsample名ファイルを利用して、inputとなる`sample.txt`を作成します。
下記のコマンドを実施してください

```
qsub pre_run.sh
```

### 5. Run the pipeline

nextflowを実装します。
`module load java/8`を先に実行しないとうまくいきません。

```
module load java/8

~/bin/nextflow run nextflow/main.nf -c run.config -resume -with-report log.01.main.html
~/bin/nextflow run nextflow/hisat2.nf -c run.config -resume -with-report log.02.hisat2.html
~/bin/nextflow run nextflow/stringtie.nf -c run.config -resume -with-report log.03.stringtie.html
~/bin/nextflow run nextflow/qc_for_R.nf -c run.config -resume -with-report log.04.qc_for_R.html
~/bin/nextflow run nextflow/summary.nf -c run.config -resume -with-report log.05.summary.html
``` 

現在1行の実施で完了するように改良中です。アップデートをお待ちください。

### 6. Output files

結果はsummaryディレクトリに作成されます。
