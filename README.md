# NF_rseq_pipeline 

## Requirements
- [Nextflow](https://www.nextflow.io/)
- [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

## prepare

### Nextflowを導入する
```
qlogin -l s_vmem=64G

以下を実⾏するとnextfileバイナリが⽣成されるので、パスを通す
curl -s https://get.nextflow.io | bash
chmod +x nextflow
```

パスは自分は~/.bash_profileに記載したが、~/.bashrc.intrに記載するのがいいのかもしれない。
```
export PATH=/home/<DOWNLOAD_PATH>/nextflow.<DOWNLOAD_VERSION>:$PATH
```

shirokaneのjavaのデフォルトはversion10。nextflowはversion 11 up to 18のjavaを使うことが推奨されているので、動かないならそこもポイントになるかもしれない。(自分は動かなかった。)

shirokaneにはjava/13が入っているので、以下を.bash_profileに追記する。
```
module load java/13
```

### Nextflow Towerのアカウントを作成
[Nextflow Tower](https://help.tower.nf/22.2/getting-started/usage/)
Sign inページから作成する。自分はGitHubアカウントとの連携で作成した。
取得したトークンは実行スクリプトである、`nf.sh`に記載する。

```
export TOWER_ACCESS_TOKEN=****
export NXF_VER=<DOWNLOAD_VERSION>
```

### ここまでのテスト

tutorial.nfを作成し、towerで監視しながら実行する。

tutorial.nf

```nextfloew
#!/usr/bin/env nextflow
params.str = 'Hello world!'
process splitLetters {
  output:
  file 'chunk_*' into letters
  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}
process convertToUpper {
  input:
  file x from letters.flatten()
  output:
  stdout result
  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}
result.view { it.trim() }
```

実行

```bash
nextflow run tutorial.nf -with-tower
```

nextflow towerのページでrunを監視できる

## run for mouse
```
git clone https://github.com/TsumaR/NF_rseq_pipeline.git
```

1. fastqファイルが入った、`fastq`という名前のフォルダをNF_rseq_pipeline直下におく。
2. `nf.sh`ファイルの下記のXXX部分を自分のものに合うように書き換える。
```
export TOWER_ACCESS_TOKEN=XXXX
export NXF_VER=XXX
```

```
qsub run.sh
```

## run for human
nextflow.configの
```
includeConfig "./config_file/mm.config"
```
を、
```
includeConfig "./config_file/hs.config"
```
に書き換える。
あとは、mouseの時と一緒。