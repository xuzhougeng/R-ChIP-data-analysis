# 如何进行R-ChIP数据分析

提高自己分析能力的一个好的方法就是重复别人文章里的分析策略，所以这里会尝试对第一篇介绍R-ChIP技术文章"R-ChIP Using Inactive RNase H Reveals Dynamic Coupling of R-loops with Transcriptional Pausing at Gene Promoters"里的所有分析进行重复，我重复所用代码会更新在我的GitHub上，地址为<https://github.com/xuzhougeng/R-ChIP-data-analysis>

选择这篇文章进行重复的理由有三点:

一：最近要探索R-loop数据分析流程

二：这篇文章的通讯作者是大牛，Xiang-Dong Fu

三：这篇文章将分析所用代码都托管在<https://github.com/Jia-Yu-Chen>

## 背景知识

高通量数据分析其中在流程上大同小异，大多数会跑流程却不知道如何分析的人，更多缺乏的是生物学基础。在阅读文献后，我整理下几个我觉得和后续数据分析有关的几个知识点:

- R-loop是一种RNA/DNA三链结构体，与基因组稳定性和转录调控有关。
- 通过电镜观察，R-loop大小在150~500bp之间。
- 硫酸氢盐测序(bisulfate sequencing)表明R-loop主要出现在基因启动子的下游。
- R-loop所在非模板链(又称编码链)具有很强的序列偏好性，计算方式为(G-C)/(G+C)

R-loop的高通量分析方法目前都是依赖于S9.6抗体捕获RNA/DNA杂合体，然后超声打断或酶切，如果后续对DNA进行测序，那就是DRIP-seq(DNA:RNA immunoprecipitation [DRIP] sequencing)，如果后续对RNA逆转成的cDNA继续测序，那就是 [DRIPc]-seq(DNA:RNA immunoprecipitation followed by cDNA conversion)。 然而酶切的分辨率不够，超声又容易破坏脆弱的R-loop结构，于是就导致目前很多文献报道有矛盾。

这篇文章就开发了一种新方法，基于RNase H的体内R-loop谱检测策略。作者构建一种没有催化活性的RNASE H1(原始的RNase H会同时靶向和基因组和线粒体基因组，因此在N端加上了NLS用于定位，在C端加上V5标签用于IP ），RNASEH1与RNA/DNA结合，超声打碎，用anti-V5抗体进行染色体免疫共沉淀(ChIP)。随后RNA/DNA杂合体转换成双链DNA(ds-DNA),  之后便是链特异性测序。

![R-loop](http://oex750gzt.bkt.clouddn.com/18-9-9/30079085.jpg)

关于链特异性测序，推荐拜读[链特异性测序那点事](http://kaopubear.top/2017-10-23-ssrnaseqbasic.html)

> “和普通的RNAseq不同，链特异性测序可以保留最初产生RNA的方向”

## 文章的分析内容

后续的分析部分会尝试完成文章中出现的分析。

1. 证明R-ChIP方法的效率和可靠度
   - 找到D210N特异而WKKD非特异的基因
   - peak的正负链比对情况
1. 序列偏好性和基因组分布
   - 比较narrow Peak和broad peak的peak大小（boxplot)
   - broad peak是否包含narrow peak（韦恩图），且一方特异的peak的信号值也弱
   - G/C skew， 统计距离summit的ATCG比例，统计背景和R-loop相比，G二聚体、G-三聚体、G四聚体..的出现比例
   - G-rich 在非模板链出现次树显著增加
   - peak在基因组上的分布(promoter proximal regions, +- Kb from TSS). 比较不同基因不同部分的差异
1. 与S9.6的R-loop结果系统性比较(基于K562 cells)
   - 比较DRIP-seq和R-ChIP的peak大小
   - 看两者的peak的overlap
   - 在JUN座位上的R-ChIP, DRIP-seq, DNase, H3K4me1, H3K4me2,H3K4me3, H3K27ac
1. 人类基因组上其他R-loop热点分析
1. R-loop诱导转录
1. R-loop的形成需要游离RNA 末端(free RNA ends)

## 准备分析环境

### 软件部分

文章中"Software and Algorithms"这部分列出了分析主要所用的软件，加上下载SRA数据所需工具和一些常用软件，一共要安装的软件如下:

- SRA Toolkit: 数据下载工具
- Bowtie2: 比对工具
- SAMtools: SAM格式处理工具
- BEDtools: BED格式处理工具
- MACS2: 比对后找peak
- R: 统计作图
- Ngsplot: 可视化工具
- Deeptools:  BAM文件分析工具, 可作图。

软件安装部分此处不介绍，毕竟如果你连软件安装都有困难，那你应该需要先学点Linux基础，或者去看[生信必修课之软件安装](https://ke.qq.com/course/310838?_bid=167&_wv=3)

### 分析项目搭建

使用`mkdir`创建项目文件夹，用于存放后续分析的所用到的数据、中间文件和结果

```bash
mkdir -p r-chip/{analysis/0-raw-data,index,scripts,results}
```

个人习惯，在项目根目录下创建了四个文件夹

- analysis: 存放原始数据、中间文件
- index: 存放比对软件索引
- scripts: 存放分析中用到的脚本
- results: 存放可用于放在文章中的结果

后续所有的操作都默认在`r-chip`下进行，除非特别说明。

### 数据下载

根据文章提供的GEO编号(GEO: GSE97072)在NCBI上检索, 按照如下步骤获取该编号下所有数据的元信息, 我将其重命名为"download_table.txt"然后上传到服务器 。

![获取数据元信息](http://oex750gzt.bkt.clouddn.com/18-9-8/88894393.jpg)

然后在`scripts`下新建一个脚本`01.sra_download.sh`

```bash
#!/bin/bash

sra_files=${1?missing input file}

sra_ids=$(egrep -o 'SRR[0-9]{5,9}' $sra_files)

mkdir -p analysis/log

for id in $sra_ids
do
    prefetch $id >> analysis/log/download.log
done
```

在项目根目录下运行`bash scritps/01.sra_download.sh ownload_table.txt`

下载的数据默认情况下存放在`~/ncbi/public/sra`, 需要用`fastq-dump`解压缩到`analysis/0-raw-data`. `fastq-dump`的使用说明见[Fastq-dump: 一个神奇的软件](https://www.jianshu.com/p/a8d70b66794c)

新建一个脚本，叫做`02.sra_uncompress.sh`，存放在`scripts`文件下，代码如下

```bash
#!/bin/bash

file_in=${1?missing input file}
SRR_IDS=$( egrep -o 'SRR[0-9]{5,9}' ${file_in})
OUT_DIR="analysis/0-raw-data"
mkdir -p ${OUT_DIR}

for id in $SRR_IDS
do
   echo "fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' $id -O ${OUT_DIR} &" | bash
done
```

然后用`bash 02.sra_uncompress.sh download_table.txt`运行。

> 如果空间不够用，可以用`egrep -o 'SRR[0-9]{5,9}' download_table.txt | xargs -i rm ~/ncbi/public/sra/{}.sra`删除SRR数据。

**注意**：这是单端测序，所以每个SRR只会解压缩出一个文件

此外还需要下载human genome (hg19)的bowtie2索引，用于后续bowtie2比对。

```bash
curl -s ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip -o index/hg19.zip &
cd index
unzip hg19.zip
```

从索引提取每条染色体的长度

```bash
bowtie2-inspect -s index/hg19 | cut -f 2,3 | grep '^chr' > index/hg19.chrom.sizes
```

## 数据分析

### 文件重命名

我们需要对下载的SRRXXXXX文件进行重命名，毕竟有意义的命名才能方便后续展示。那么，应该如何做呢？

首先，你需要将[GSE97072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97072)页面的中Samples这部分的内容复制到一个文本文件中（我将其命名为sample_name.txt)，分为两列，第一列是GSM编号，第二列是样本的命名。

![sample name](http://oex750gzt.bkt.clouddn.com/18-9-15/73109617.jpg)

> 注：这里面有一个希腊字符在不同系统表示有所不同，所以我在复制之后手动修改。

此外，你还需要将里面的`(-)`和`(+)`进行替换，因为括号在shell里有特殊含义， 为了保证命名的连贯性，我将`(-)`替换成`-neg`，将`(+)`替换成`-pos`.

```bash
 sed -i 's/(-)/-neg/; s/(+)/-pos/' sample_name.txt
```

随后，我们需要从download_table.txt中提取出SRR编号和GSM编号的对应的关系，这个需要用到Linux的文本处理命令

```bash
paste <(egrep -o 'GSM[0-9]{6,9}' download_table.txt ) <(egrep -o 'SRR[0-9]{6,9}' download_table.txt) >  gsm_srr.txt
join <(sort gsm_srr.txt ) <(head -n 24 sample_name.txt | sort) > gsm_srr_sample_name.txt
```

> 注: 样本一共有25行，而我们只下载了24个数据，所以删除了最后一行的"K562-WKKD-V5chip"

最后，就是根据生成的文件对样本进行重命名了。

```bash
awk '{print "rename " $2 " " $3   " analysis/0-raw-data/" $2 "*.gz"}' gsm_srr_sample_name.txt |bash
```

这行代码看起来有点复杂，但是其实做的事情就是构建了一系列`rename`的命令行，然后在bash下执行。

### R-ChIP分析

#### 数据预处理

这一部分主要是将序列比对后原来的参考基因组上，标记重复，并且去掉不符合要求的比对。

让我们先写一个比对脚本将序列比对到参考基因组上，脚本命名为`03.r_chip_align.sh`，存放在`scripts`下

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

# configuration
threads=8
index=index/hg19
FQ_DIR="analysis/0-raw-data"
ALIGN_DIR="analysis/2-read-align"
LOG_DIR="analysis/log"
TMP_DIR="analysis/tmp"

mkdir -p ${ALIGN_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${TMP_DIR}

samples=${1?missing sample file}

exec 0< $samples
# alignment
while read id;
do
    if [ ! -f ${FQ_DIR}/$id.bt2.done ]
    then
    echo "bowtie2 --very-sensitive-local --mm -p $threads -x $index -U ${FQ_DIR}/$id.fastq.gz 2> ${LOG_DIR}/$id.bt2.log | \
    samtools sort -@ 2 -m 1G -T ${TMP_DIR}/${id} -o ${ALIGN_DIR}/${id}.sort.bam \
    && touch ${ALIGN_DIR}/$id.bt2.done" | bash
    fi
done
```

这个脚本的前半部分都在定义各种变量，而`#alignment`标记的后半部分则是实际运行的比对命令

> 你会发现我的实际运行部分脚本有点奇怪，我没有直接运行比对，而是用`echo`通过管道传递给了bash执行。
>
> 这样做的原因是, 当我用`sed -i 's/| bash/#| bash/ 03.r_chip_align.sh'`将`| bash`替换成`#| bash`后，运行这个脚本就可以检查代码是否正确，然后再用`sed-i s/#| bash/| bash/03.r_chip_align.sh`修改回来 。 后面的脚本同理。

接着再写一个脚本用于标记重复`04.markdup.sh`，也是放在`scripts`下

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

# configuration
threads=8
index=index/hg19
FQ_DIR="analysis/0-raw-data"
OUT_DIR="analysis/2-read-align"
LOG_DIR="analysis/log"
TMP_DIR="analysis/tmp"

mkdir -p ${OUT_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${TMP_DIR}

samples=${1?missing sample file}

exec 0< $samples
# mark duplication
while read id;
do
    if [ ! -f ${OUT_DIR}/${id}.mkdup.done ]
    then
    echo "sambamba markdup -t $threads ${OUT_DIR}/${id}.sort.bam ${OUT_DIR}/${id}.mkdup.bam \
    && touch ${OUT_DIR}/${id}.mkdup.done" | bash
    fi
done
```

再准备一个脚本用于去除不符合要求的比对，命名为`05.bam_filter.sh`

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

# configuration
threads=8
index=index/hg19
FQ_DIR="analysis/0-raw-data"
ALIGN_DIR="analysis/2-read-align"
LOG_DIR="analysis/log"
TMP_DIR="analysis/tmp"

mkdir -p ${ALIGN_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${TMP_DIR}

samples=${1?missing sample file}

exec 0< $samples
# filter
while read id;
do
    if [ ! -f ${ALIGN_DIR}/${id}.flt.done ]
    then
    echo "
    samtools view -@ threads -bF 1804 -q 30 ${ALIGN_DIR}/${id}.mkdup.bam -o ${ALIGN_DIR}/${id}.flt.bam\
    && samtools index ${ALIGN_DIR}/${id}.flt.bam \
    && touch ${ALIGN_DIR}/${id}.flt.done "| bash
    fi
done
```

最后新建一个`samples01.txt`用于存放将要处理的样本名

```bash
HKE293-D210N-V5ChIP-Rep1
HKE293-D210N-Input-Rep1
HKE293-D210N-V5ChIP-Rep2
HKE293-D210N-Input-Rep2
HKE293-D210N-V5ChIP-Rep3
HKE293-D210N-Input-Rep3
HKE293-WKKD-V5ChIP
HKE293-WKKD-Input
HKE293-delta-HC-V5ChIP
```

这样子就可以依次运行脚本了

- 比对: `bash scripts/03.r_chip_align.sh samples01.txt`
- 标记重复：`bash scripts/04.bam_markdup.sh samples01.txt`
- 去除不符合要求的比对: `bash scripts/05.bam_filter.sh samples01.txt`

处理完之后可以对每个样本都进行一次统计，包括如下信息：

- 处理前的原始reads数
- 处理后对唯一比对reads数
- 唯一比对reads数所占原始reads数的比例

这个工作同样可以用shell脚本完成, 脚本为`06.sample_align_stat.sh`

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?missing sample file}
threads=8
OUT_DIR="analysis/2-read-align"

echo -e "Experiment \t Raw Reads \t Uniquely mapped Reads \t ratio"
exec 0< $samples

while read sample;
do
     total=$( samtools view -@ ${threads} -c ${ALIGN_DIR}/${sample}.sort.bam )
     unique=$( samtools view -@ ${threads} -c ${ALIGN_DIR}/${sample}.flt.bam )
     ratio=$( echo "scale=2; 100 * $unique / $total " | bc )
     echo -e "$sample \t $total \t $unique \t $ratio %"
done
```

运行方法是`bash scripts/06.sample_align_stat.sh samples01.txt > results/library_stat.txt`，运行结果如下，和原本的Table S2对比，你会发现结果基本一致，有出入的地方我推测是标记重复这一步所用软件不同。

| Experiment               | Raw   Reads | Uniquely mapped Reads | ratio  |
| ------------------------ | ----------- | --------------------- | ------ |
| HKE293-D210N-V5ChIP-Rep1 | 22405416    | 6443979               | 28.76% |
| HKE293-D210N-Input-Rep1  | 60302237    | 25673307              | 42.57% |
| HKE293-D210N-V5ChIP-Rep2 | 17763614    | 11778533              | 66.30% |
| HKE293-D210N-Input-Rep2  | 11131443    | 8553097               | 76.83% |
| HKE293-D210N-V5ChIP-Rep3 | 8799855     | 5640375               | 64.09% |
| HKE293-D210N-Input-Rep3  | 4529910     | 3209275               | 70.84% |
| HKE293-WKKD-V5ChIP       | 12734577    | 8612940               | 67.63% |
| HKE293-WKKD-Input        | 8830478     | 6643507               | 75.23% |
| HKE293-delta-HC-V5ChIP   | 25174573    | 9252009               | 36.75% |

#### BAM相关性评估

上一步得到各个样本的BAM文件之后，就可以在全基因组范围上看看这几个样本之间是否有差异。也就是先将基因组分成N个区间，然后用统计每个区间上比对上的read数。

脚本`scripts/07.genome_bin_read_coverage.sh`如下

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?missing sample file}
chromsize=${2:-index/hg19.chrom.sizes}
size=${3:-3000}
ALIGN_DIR="analysis/2-read-align"
COV_DIR="analysis/3-genome-coverage"

mkdir -p ${COV_DIR}

exec 0< $samples

while read sample
do
    bedtools makewindows -g $chromsize -w $size | \
      bedtools intersect -b ${ALIGN_DIR}/${sample}.flt.bam -a - -c -bed > ${COV_DIR}/${sample}.ReadsCoverage
done

```

准备一个输入文件存放待处理样本的前缀，然后运行脚本`bash scripts/07.genome_bin_read_coverage.sh samples_rep.txt`

```bash
HKE293-D210N-Input-Rep1
HKE293-D210N-Input-Rep2
HKE293-D210N-Input-Rep3
HKE293-D210N-V5ChIP-Rep1
HKE293-D210N-V5ChIP-Rep2
HKE293-D210N-V5ChIP-Rep3
```

最后将得到的文件导入到R语言中进行作图，使用的是基础绘图系统的光滑散点图(smoothScatter)。

```r
files_list <- list.files("r-chip/analysis/3-genome-coverage","ReadsCoverage")
files_path <- file.path("r-chip/analysis/3-genome-coverage",files_list)

input_rep1 <- read.table(files_path[1], sep='\t')
input_rep2 <- read.table(files_path[2], sep='\t')
input_rep3 <- read.table(files_path[3], sep='\t')
chip_rep1 <- read.table(files_path[4], sep='\t')
chip_rep2 <- read.table(files_path[5], sep='\t')
chip_rep3 <- read.table(files_path[6], sep='\t')

pw_plot <- function(x, y, 
                    xlab="x",
                    ylab="y", ...){
  log2x <- log2(x)
  log2y <- log2(y)
  smoothScatter(log2x,log2y,
                cex=1.2,
                xlim=c(0,12),ylim=c(0,12),
                xlab=xlab,
                ylab=ylab)
  text(3,10,paste("R = ",round(cor(x,y),2),sep=""))
}

par(mfrow=c(2,3))
pw_plot(chip_rep1[,4], chip_rep2[,4],
        xlab = "Rep 1 (Log2 Tag Counts)",
        ylab = "Rep 2 (Log2 Tag Counts)")

pw_plot(chip_rep1[,4], chip_rep3[,4],
        xlab = "Rep 1 (Log2 Tag Counts)",
        ylab = "Rep 2 (Log2 Tag Counts)")

pw_plot(chip_rep2[,4], chip_rep3[,4],
        xlab = "Rep 1 (Log2 Tag Counts)",
        ylab = "Rep 2 (Log2 Tag Counts)")

pw_plot(input_rep1[,4], input_rep2[,4],
        xlab = "Rep 1 (Log2 Tag Counts)",
        ylab = "Rep 2 (Log2 Tag Counts)")

pw_plot(input_rep1[,4], input_rep3[,4],
        xlab = "Rep 1 (Log2 Tag Counts)",
        ylab = "Rep 3 (Log2 Tag Counts)")

pw_plot(input_rep2[,4], input_rep3[,4],
        xlab = "Rep 2 (Log2 Tag Counts)",
        ylab = "Rep 3 (Log2 Tag Counts)")
```

下图上为D210的R-ChIP三个重复间的相关性，下为Input三个重复间的相关性

![correlationship](http://oex750gzt.bkt.clouddn.com/18-9-17/30639991.jpg)

> 这种多个BAM文件之间相关性衡量，其实也可以用`deepTools`的`plotCorrelation`画出来，但是我觉得应该没有R语言画 的好看。

**看图说话**： 由于同一个样本间的BAM文件具有很强的相关性，因此可以将这些样本合并起来，这样子在基因组浏览器上就可以只用一个轨（tracks）

#### BigWig可视化

虽然可以直接用BAM也行对比对结果进行可视化展示，但是一般BAM文件文件太大，不方便传输，所以需要转换成bigwig格式，这样子在基因组浏览器（例如IGV, UCSC Browser, JBrowse）上方便展示。

如下的代码的目的就是先合并BAM，然后转换成BigWig，拆分成正链和反链进行保存

```bash
#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?please provied sample file}
threads=${2-8}
bs=${3-50}

ALIGN_DIR="analysis/2-read-align"
BW_DIR="analysis/4-normliazed-bw"

mkdir -p $BW_DIR

exec 0< $samples
cd $ALIGN_DIR

while read sample
do
    bams=$(ls ${sample}*flt.bam | tr '\n' ' ')
    if [ ! -f ../$(basename ${BW_DIR})/${sample}.tmp.bam ]; then
    echo "samtools merge -f -@ ${threads} ../$(basename ${BW_DIR})/${sample}.tmp.bam $bams  &&
          samtools index ../$(basename ${BW_DIR})/${sample}.tmp.bam" | bash
    fi
done

cd ../../
exec 0< $samples

cd ${BW_DIR}
while read sample
do
    bamCoverage -b ${sample}.tmp.bam -o ${sample}_fwd.bw -of bigwig  \
      --filterRNAstrand forward --binSize ${bs} --normalizeUsing CPM --effectiveGenomeSize 2864785220 \
      --extendReads 150 -p ${threads} 2> ../log/${sample}_fwd.log
    bamCoverage -b ${sample}.tmp.bam -o ${sample}_rvs.bw -of bigwig  \
      --filterRNAstrand reverse --binSize ${bs} --normalizeUsing CPM --effectiveGenomeSize 2864785220 \
      --extendReads 150 -p ${threads} 2> ../log/${sample}_rvs.log
    rm -f ${sample}.tmp.bam ${sample}.tmp.bam.bai
done

```

得到的BW文件你可以在IGV上初步看看，比如说检查下文章Figure 1(E)提到的CIRH1A基因

![CIRH1A](http://oex750gzt.bkt.clouddn.com/18-9-18/90750745.jpg)

发表文章时肯定不能用上图，我们可以用R的`Gviz`进行展示下

```r
library(Gviz)

#下面将scale等track写入tracklist
tracklist<-list()
itrack <- IdeogramTrack(genome = "hg19", chromosome = 'chr16',outline=T)
tracklist[["itrack"]]<-itrack

# 读取BigWig
bw_file_path <- "C:/Users/DELL/Desktop/FigureYa/R-ChIP/"
bw_file_names <- list.files(bw_file_path, "*.bw")
bw_files <- file.path("C:/Users/DELL/Desktop/FigureYa/R-ChIP/",
                      bw_file_names)

tracklist[['D210-fwd']] <- DataTrack(range = bw_files[3],
                              genome="hg19",
                              type="histogram",
                              name='D210 + ',
                              ylim=c(0,4),
                              col.histogram="#2167a4",
                              fill.histogram="#2167a4")
tracklist[['D210-rvs']] <- DataTrack(range = bw_files[4],
                                     genome="hg19",
                                     type="histogram",
                                     name='D210 - ',
                                     ylim=c(4,0),
                                     col.histogram="#eb1700",
                                     fill.histogram="#eb1700")

tracklist[['WKDD-fwd']] <- DataTrack(range = bw_files[9],
                                     genome="hg19",
                                     type="histogram",
                                     name='WKDD + ',
                                     ylim=c(0,4),
                                     ylab=2,
                                     col.histogram="#2167a4",
                                     fill.histogram="#2167a4")
tracklist[['WKDD-rvs']] <- DataTrack(range = bw_files[10],
                                     genome="hg19",
                                     type="histogram",
                                     name=' WKDD -',
                                     ylim=c(4,0),
                                     col.histogram="#eb1700",
                                     fill.histogram="#eb1700",
                                     showAxis=TRUE)

#写入比例尺
scalebar <- GenomeAxisTrack(scale=0.25,
                            col="black",
                            fontcolor="black",
                            name="Scale",
                            labelPos="above",showTitle=F)
tracklist[["scalebar"]]<-scalebar

# 画图
plotTracks(trackList = tracklist,
           chromosome = "chr16",
           from =  69141913, to= 69205033,
           background.panel = "#f6f6f6",
           background.title = "white",
           col.title="black",col.axis="black",
           rot.title=0,cex.title=0.5,margin=10,title.width=1.75,
           cex.axis=1
           )
```

![原文Fig1E](http://oex750gzt.bkt.clouddn.com/18-9-18/62177187.jpg)

后续用AI修改下坐标轴，就几乎和原图差不多了。此外可能还要调整之前脚本的bin size， 使得整体更加平滑。

**看图说话**： 由于WKKD在D210N的基础上继续突变了几个位点，使得原本只丧失了RNASEH1的催化活性的D210N进一步丧失了结合到RNA/DNA杂合体的能力，在实验中就可以作为 **负对照**。上图就是其中一个有代表性的区间。

#### Peak Calling

关于MACS2的使用方法， 我写了[如何使用MACS进行peak calling](https://www.jianshu.com/p/6a975f0ea65a)详细地介绍了它的参数。按照默认参数分别找narrow peak 和 broad peak

```bash
cd analysis
# HKE293 D210N
macs2 callpeak -g hs --keep-dup all -f BAM -t 2-read-align/HKE293-D210N-V5ChIP-Rep*.flt.bam -c 2-read-align/HKE293-D210N-Input-Rep*.flt.bam --outdir 5-peak-calling/narrow -n HKE293-D210 2> log/HKE293-D210.macs2.narrow.log &
macs2 callpeak -g hs --keep-dup all --broad -f BAM -t 2-read-align/HKE293-D210N-V5ChIP-Rep*.flt.bam -c 2-read-align/HKE293-D210N-Input-Rep*.flt.bam --outdir 5-peak-calling/broad -n HKE293-D210 2> log/HKE293-D210.macs2.broad.log &
# HEK293 WKKD
macs2 callpeak -g hs --keep-dup all -f BAM -t 2-read-align/HKE293-WKKD-V5ChIP.flt.bam -c 2-read-align/HKE293-WKKD-Input.flt.bam --outdir 5-peak-calling/narrow/ -n HKE293-WKKD 2> log/HKE293-WKKD.macs2.narrow.log &
macs2 callpeak -g hs --keep-dup all --broad -f BAM -t 2-read-align/HKE293-WKKD-V5ChIP.flt.bam -c 2-read-align/HKE293-WKKD-Input.flt.bam --outdir 5-peak-calling/broad/ -n HKE293-WKKD 2> log/HKE293-WKKD.macs2.broad.log &
```

文章中对找到peak进行了一次筛选，标准是大于5倍富集和q-value 小于或等于0.001（broad peak则是0.0001），最后文章写着在D210N有12,906peak，然后剔除WKKD里的peak，还有12,507个。

我找到的HKE293-D210的原始narrowPeak数为15,026, 按照作者的标准筛选后只剩下6639，发现负对照的WKKD只有2个peak，从D210过滤后基础上又剔除了一个。对于HKE2930-D210的原始broadPeak数为39278，过滤之后只剩2358了。

这似乎表明在我的分析流程下，原标准有点过于严格了，因此我这里使用3倍。

```bash
cd analysis/5-peak-calling/
fc=3
# narrow peak
awk -v fc=$fc '$7 >= fc && $9 >=3' narrow/HKE293-D210_peaks.narrowPeak > narrow/HKE293-D210_peaks.narrowPeak.tmp
bedtools subtract -A -a narrow/HKE293-D210_peaks.narrowPeak.tmp -b narrow/HKE293-WKKD_peaks.narrowPeak > narrow/HKE293-D210_flt.narrowPeak
rm -f narrow/HKE293-D210_peaks.narrowPeak.tmp
# broad peak
awk -v fc=$fc '$7 >= fc && $9 >=4' broad/HKE293-D210_peaks.broadPeak > broad/HKE293-D210_peaks.broadPeak.tmp
bedtools subtract -A -a broad/HKE293-D210_peaks.broadPeak.tmp -b broad/HKE293-WKKD_peaks.broadPeak > broad/HKE293-D210_flt.broadPeak
rm -f broad/HKE293-D210_peaks.broadPeak.tmp
```

之后可以用过滤后的D210N的narrowPeak和broadPeak的peak长度进行描述性统计分析，然后用箱线图展示其大小分布。

```r
library(data.table)
library(ggplot2)

narrowPeak <- fread(file="HKE293-D210_flt.narrowPeak",
                   sep="\t", header = F)
broadPeak <- fread(file="HKE293-D210_flt.broadPeak",
                    sep="\t", header = F)
peak_size <- log10(c(narrowPeak$V3 - narrowPeak$V2, broadPeak$V3 - broadPeak$V2 ))

peak_from <- factor(rep(c('Narrow','Broad'), times=c(nrow(narrowPeak),nrow(broadPeak)) ),
                    levels=c("Narrow","Broad"))

peak_df <- data.frame(size=peak_size, from=peak_from)

my_clear_theme = theme_bw() + 
  theme(panel.grid.major = element_line(colour="NA"),
        panel.grid.minor = element_line(colour="NA"),
        panel.background = element_rect(fill="NA"),
        panel.border = element_rect(colour="black", fill=NA),
        legend.background = element_blank())

ggplot(peak_df,aes(x=from,y=size,col=from)) + 
  geom_boxplot(notch = T,outlier.size = .5) +
  scale_y_continuous(breaks=c(log10(100),log10(300),log10(400),log10(450),log10(500),log10(1000),log10(2000)),
                     labels = c(100,300,400,450,500,1000,2000)) +
  coord_cartesian(ylim=c(2,3.5)) + labs(col="Peak Size") +
  my_clear_theme + xlab("") + ylab("Peak Size (bp)")
ggsave("Fig2A_boxplot.pdf", width=4,height=6,units = "in")
```

![Fig2A](http://oex750gzt.bkt.clouddn.com/18-9-19/37564618.jpg)

原文说自己的narrow peak 长度的中位数是199bp, broad peak的长度中位数是318 bp，和电镜观察的150–500 bp一致。然而我过滤后的peak的中位数，除了narrow peak的中位数是380 bp勉强在150-500之间，broad peak的中位数已经快突破天际了。

> 我觉得这或许和MACS2的版本有点关系，也可能是作者其实是分别用实验组和对照组找peak，然后进行把peak进行合并？也有可能是作者把BAM转换成了BED然后分析。这就是留给大家去验证，但是narrow peak的中位数和电镜结果一致就证明了这个方法是比较成功啦

还可以对Broad peak 和Narrow peak进行比较，看看有多少共同peak和特异性peak。

```bash
#!/bin/bash
#!/bin/bash

peak1=${1?first peak}
peak2=${2?second peak}
outdir=${3?output dir}

mkdir -p ${outdir}

bedtools intersect -a $peak1 -b $peak2 > ${outdir}/$(basename ${peak1%%.*}).common.peak
common=$(wc -l ${outdir}/$(basename ${peak1%%.*}).common.peak)
echo "$common"

bedtools subtract -A -a $peak1 -b $peak2 > ${outdir}/$(basename ${peak1%%.*}).only.peak
left=$(wc -l ${outdir}/$(basename ${peak1%%.*}).only.peak)

echo "$left"

bedtools subtract -A -b $peak1 -a $peak2 > ${outdir}/$(basename ${peak2%%.*}).only.peak
right=$(wc -l ${outdir}/$(basename ${peak2%%.*}).only.peak)

echo "$right"
```

运行:`bash scripts/09.peak_compare.sh analysis/5-peak-calling/narrow/HKE293-D210_flt.narrowPeak analysis/5-peak-calling/broad/HKE293-D210_flt.broadPeak analysis/5-peak-calling/peak_compare > results/narrow_vs_broad.txt`

然后得到的数值就可以丢到R语言画图了，这个结果和Figure S2A一致。

```r
library(VennDiagram)
grid.newpage()
venn.plot <- VennDiagram::draw.pairwise.venn(area1=6401 + 7165,
                                area2=286 + 7165,
                                cross.area = 7165,
                                category = c("Narrow specific","Broad specific"),
                                scaled = F,
                                fill=c("#c7d0e7","#f4babf"),
                                lty='blank',
                                cat.pos = c(180,0),#360度划分
                                rotation.degree = 90 #整体旋转90
                                )
grid.draw(venn.plot)

```

![venn plot](http://oex750gzt.bkt.clouddn.com/18-9-20/29733091.jpg)

我们可以对这三类peak对应区间内的信号进行一下比较。统计信号强度的工具是`bigwigSummary`，来自于[ucscGenomeBrowser](https://github.com/ucscGenomeBrowser/kent)工具集。

脚本为`scripts/10.peak_signal_summary.sh`

```bash
#!/bin/bash

bed=${1?bed file}
fwd=${2?forwad bigwig}
rvs=${3?reverse bigwig}

exec 0< $bed

while read region
do
    chr=$( echo $region | cut -d ' ' -f 1)
    sta=$( echo $region | cut -d ' ' -f 2)
    end=$( echo $region | cut -d ' ' -f 3)
    (bigWigSummary $fwd $chr $sta $end 1 ; bigWigSummary $rvs $chr $sta $end 1 )| awk -v out=0 '{out=out+$1} END{print out/2}'
done
```

分别统计三类peak的信号的信号

```bash
bash scripts/10.peak_signal_summary.sh analysis/5-peak-calling/peak_compare/HKE293-D210_flt.broadPeak.only.bed analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_fwd.bw analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_rvs.bw > results/HKE293-D210N-V5ChIP.broadSignal
bash scripts/10.peak_signal_summary.sh analysis/5-peak-calling/peak_compare/HKE293-D210_flt.narrowPeak.only.bed analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_fwd.bw analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_rvs.bw > results/HKE293-D210N-V5ChIP.narrowSignal
bash scripts/10.peak_signal_summary.sh analysis/5-peak-calling/peak_compare/HKE293-D210_flt.common.peak analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_fwd.bw analysis/4-normliazed-bw/HKE293-D210N-V5ChIP_rvs.bw > results/HKE293-D210N-V5ChIP.commonSignal
```

在R语言中绘图展示

```bash
# peak signal
commonSignal <- fread("./HKE293-D210N-V5ChIP.commonSignal")
narrowSignal <- fread("./HKE293-D210N-V5ChIP.narrowSignal")
broadSignal <-  fread("./HKE293-D210N-V5ChIP.broadSignal")

data <- data.frame(signal=c(commonSignal$V1, narrowSignal$V1, broadSignal$V1),
                   type=factor(rep(c("common","narrow","broad"),
                            times=c(nrow(commonSignal),nrow(narrowSignal), nrow(broadSignal)))) )

wilcox.test(data[data$type=="common",1],data[data$type=="broad",1])
wilcox.test(data[data$type=="common",1],data[data$type=="narrow",1])


p <- ggplot(data,aes(x=type,y=signal,col=type)) + 
  geom_boxplot(outlier.colour = "NA") +
  coord_cartesian(ylim=c(0,3.5)) +
  theme(panel.grid.major = element_line(colour="NA"),
        panel.grid.minor = element_line(colour="NA"),
        panel.background = element_rect(fill="NA"),
        panel.border = element_rect(colour="black", fill=NA)) +
  theme(axis.text.x = element_blank()) +
  xlab("") + ylab("R-loop Signal") + labs(col="") 


p1 <- p + annotate("segment", x=1,xend=1,y=0.75,yend=2,size=0.5) + 
  annotate("segment", x=2,xend=2,y=1.25,yend=2,size=0.5) +
  annotate("segment", x=1,xend=2,y=2,yend=2,size=0.5) +
  annotate("text", x= 1.5, y=2.1, label="***") 

p2 <- p1 + annotate("segment", x=3,xend=3,y=0.65,yend=3,size=0.5) + 
  annotate("segment", x=2,xend=2,y=2.2,yend=3,size=0.5) +
  annotate("segment", x=2,xend=3,y=3,yend=3,size=0.5) +
  annotate("text", x= 2.5, y=3.2, label="***") 

p2
```

![signal comparison](http://oex750gzt.bkt.clouddn.com/18-9-20/36368977.jpg)

**看图说话**： 无论是narrow peak 特异区间的信号，还是broad peak 特异区间的信号都显著性低于其共有区间的信号。

