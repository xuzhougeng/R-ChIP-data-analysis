如何进行R-ChIP数据分析

这应该是全网能搜索到的第一篇介绍R-ChIP数据分析的中文教程，本篇教程所用代码更新在我的GitHub上，地址为https://github.com/xuzhougeng/R-ChIP-data-analysis

教程会对第一篇介绍R-ChIP技术文章"R-ChIP Using Inactive RNase H Reveals Dynamic Coupling of R-loops with Transcriptional Pausing at Gene Promoters"里的所有分析进行重复。

选择这篇文章进行重复的理由有三点:

一：最近要探索R-loop数据分析流程

二：这篇文章的通讯作者是大牛，Xiang-Dong Fu

三：这篇文章将分析所用代码都托管在https://github.com/Jia-Yu-Chen

背景知识

高通量数据分析其中在流程上大同小异，大多数会跑流程却不知道如何分析的人，更多缺乏的是生物学基础。在阅读文献后，我整理下几个我觉得和后续数据分析有关的几个知识点:

- R-loop是一种RNA/DNA三链结构体，与基因组稳定性和转录调控有关。
- 通过电镜观察，R-loop大小在150~500bp之间。
- 硫酸氢盐测序(bisulfate sequencing)表明R-loop主要出现在基因启动子的下游。
- R-loop所在非模板链(又称编码链)具有很强的序列偏好性，计算方式为(G-C)/(G+C)

R-loop的高通量分析方法目前都是依赖于S9.6抗体捕获RNA/DNA杂合体，然后超声打断或酶切，如果后续对DNA进行测序，那就是DRIP-seq(DNA:RNA immunoprecipitation [DRIP] sequencing)，如果后续对RNA逆转成的cDNA继续测序，那就是 [DRIPc]-seq(DNA:RNA immunoprecipitation followed by cDNA conversion)。 然而酶切的分辨率不够，超声又容易破坏脆弱的R-loop结构，于是就导致目前很多文献报道有矛盾。

这篇文章就开发了一种新方法，基于RNase H的体内R-loop谱检测策略。作者构建一种没有催化活性的RNASE H1(原始的RNase H会同时靶向和基因组和线粒体基因组，因此在N端加上了NLS用于定位，在C端加上V5标签用于IP ），RNASEH1与RNA/DNA结合，超声打碎，用anti-V5抗体进行染色体免疫共沉淀(ChIP)。随后RNA/DNA杂合体转换成双链DNA(ds-DNA),  之后便是链特异性测序。



关于链特异性测序，推荐拜读链特异性测序那点事

“和普通的RNAseq不同，链特异性测序可以保留最初产生RNA的方向”

准备分析环境

软件部分

文章中"Software and Algorithms"这部分列出了分析主要所用的软件，加上下载SRA数据所需工具和一些常用软件，一共要安装的软件如下:

- SRA Toolkit: 数据下载工具
- Bowtie2: 比对工具
- SAMtools: SAM格式处理工具
- BEDtools: BED格式处理工具
- MACS2: 比对后找peak
- R: 统计作图
- Ngsplot: 可视化工具
- Deeptools:  BAM文件分析工具, 可作图。

软件安装部分此处不介绍，毕竟如果你连软件安装都有困难，那你应该需要先学点Linux基础，或者去看生信必修课之软件安装

分析项目搭建

使用mkdir创建项目文件夹，用于存放后续分析的所用到的数据、中间文件和结果

    mkdir -p r-chip/{analysis/0-raw-data,index,scripts,results}

个人习惯，在项目根目录下创建了四个文件夹

- analysis: 存放原始数据、中间文件
- index: 存放比对软件索引
- scripts: 存放分析中用到的脚本
- results: 存放可用于放在文章中的结果

后续所有的操作都默认在r-chip下进行，除非特别说明。

数据下载

根据文章提供的GEO编号(GEO: GSE97072)在NCBI上检索, 按照如下步骤获取该编号下所有数据的元信息, 我将其重命名为"download_table.txt"然后上传到服务器 。



然后在scripts下新建一个脚本01.sra_download.sh

    #!/bin/bash
    
    sra_files=$1
    
    sra_ids=$(egrep -o 'SRR[0-9]{5,9}' $sra_files)
    
    mkdir -p analysis/log
    
    for id in $sra_ids
    do
        prefetch $id >> analysis/log/download.log
    done

在项目根目录下运行bash scritps/01.sra_download.sh ownload_table.txt

下载的数据默认情况下存放在~/ncbi/public/sra, 需要用fastq-dump解压缩到analysis/0-raw-data. fastq-dump的使用说明见Fastq-dump: 一个神奇的软件

新建一个脚本，叫做02.sra_uncompress.sh，存放在scripts文件下，代码如下

    #!/bin/bash
    
    file_in=$1
    SRR_IDS=$( egrep -o 'SRR[0-9]{5,9}' ${file_in})
    OUT_DIR="analysis/0-raw-data"
    mkdir -p ${OUT_DIR}
    
    for id in $SRR_IDS
    do
       echo "fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' $id -O ${OUT_DIR} &" | bash
    done

然后用bash 02.sra_uncompress.sh download_table.txt运行。

如果空间不够用，可以用egrep -o 'SRR[0-9]{5,9}' download_table.txt | xargs -i rm ~/ncbi/public/sra/{}.sra删除SRR数据。

注意：这是单端测序，所以每个SRR只会解压缩出一个文件

此外还需要下载human genome (hg19)的bowtie2索引，用于后续bowtie2比对。

    curl -s ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip -o index/hg19.zip &
    cd index
    unzip hg19.zip

数据分析

文件重命名

我们需要对下载的SRRXXXXX文件进行重命名，毕竟有意义的命名才能方便后续展示。那么，应该如何做呢？

首先，你需要将GSE97072页面的中Samples这部分的内容复制到一个文本文件中（我将其命名为sample_name.txt)，分为两列，第一列是GSM编号，第二列是样本的命名。



注：这里面有一个希腊字符在不同系统表示有所不同，所以我在复制之后手动修改。

此外，你还需要将里面的(-)和(+)进行替换，因为括号在shell里有特殊含义， 为了保证命名的连贯性，我将(-)替换成-neg，将(+)替换成-pos.

     sed -i 's/(-)/-neg/; s/(+)/-pos/' sample_name.txt

随后，我们需要从download_table.txt中提取出SRR编号和GSM编号的对应的关系，这个需要用到Linux的文本处理命令

    paste <(egrep -o 'GSM[0-9]{6,9}' download_table.txt ) <(egrep -o 'SRR[0-9]{6,9}' download_table.txt) >  gsm_srr.txt
    join <(sort gsm_srr.txt ) <(head -n 24 sample_name.txt | sort) > gsm_srr_sample_name.txt

注: 样本一共有25行，而我们只下载了24个数据，所以删除了最后一行的"K562-WKKD-V5chip"

最后，就是根据生成的文件对样本进行重命名了。

    awk '{print "rename " $2 " " $3   " analysis/0-raw-data/" $2 "*.gz"}' gsm_srr_sample_name.txt |bash

这行代码看起来有点复杂，但是其实做的事情就是构建了一系列rename的命令行，然后在bash下执行。

R-ChIP分析

数据预处理

这一部分主要是将序列比对后原来的参考基因组上，标记重复，并且去掉不符合要求的比对。

让我们先写一个比对脚本将序列比对到参考基因组上，脚本命名为03.r_chip_align.sh，存放在scripts下

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
    
    samples=$1
    
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

这个脚本的前半部分都在定义各种变量，而#alignment标记的后半部分则是实际运行的比对命令

你会发现我的实际运行部分脚本有点奇怪，我没有直接运行比对，而是用echo通过管道传递给了bash执行。

这样做的原因是, 当我用sed -i 's/| bash/#| bash/ 03.r_chip_align.sh'将| bash替换成#| bash后，运行这个脚本就可以检查代码是否正确，然后再用sed-i s/#| bash/| bash/03.r_chip_align.sh修改回来 。 后面的脚本同理。

接着再写一个脚本用于标记重复04.markdup.sh，也是放在scripts下

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
    
    samples=$1
    
    exec 0< $samples
    # mark duplication
    while read id;
    do
        if [ ! -f ${ALIGN_DIR}/${id}.mkdup.done ]
        then
        echo "sambamba markdup -t $threads ${ALIGN_DIR}/${id}.sort.bam ${ALIGN_DIR}/${id}.mkdup.bam \
        && touch ${ALIGN_DIR}/${id}.mkdup.done" | bash
        fi
    done

再准备一个脚本用于去除不符合要求的比对，命名为05.bam_filter.sh

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
    
    samples=$1
    
    exec 0< $samples
    # filter
    while read id;
        if [ ! -f ${ALIGN_DIR}/${id}.flt.done ]
        then
        echo "
        samtools view -@ threads -bF 1804 -q 30 ${ALIGN_DIR}/${id}.mkdup.bam -o ${ALIGN_DIR}/${id}.flt.bam\
        && samtools index ${ALIGN_DIR}/${id}.flt.bam \
        && touch ${ALIGN_DIR}/${id}.flt.done "| bash
        fi
    done

最后新建一个samples01.txt用于存放将要处理的样本名

    HKE293-D210N-V5ChIP-Rep1
    HKE293-D210N-Input-Rep1
    HKE293-D210N-V5ChIP-Rep2
    HKE293-D210N-Input-Rep2
    HKE293-D210N-V5ChIP-Rep3
    HKE293-D210N-Input-Rep3
    HKE293-WKKD-V5ChIP
    HKE293-WKKD-Input
    HKE293-delta-HC-V5ChIP

这样子就可以以此运行脚本了

- 比对: bash scripts/03.r_chip_align.sh samples01.txt
- 标记重复：bash scripts/04.bam_markdup.sh samples01.txt
- 去除不符合要求的比对: bash scripts/05.bam_filter.sh samples01.txt



下游分析任务

1. 证明R-ChIP方法的效率和可靠度
   - 找到D210N特异而WKKD非特异的基因
   - peak的正负链比对情况
2. 序列偏好性和基因组分布
   - 比较narrow Peak和broad peak的peak大小（boxplot)
   - broad peak是否包含narrow peak（韦恩图），且一方特异的peak的信号值也弱
   - G/C skew， 统计距离summit的ATCG比例，统计背景和R-loop相比，G二聚体、G-三聚体、G四聚体..的出现比例
   - G-rich 在非模板链出现次树显著增加
   - peak在基因组上的分布(promoter proximal regions, +- Kb from TSS). 比较不同基因不同部分的差异
3. 与S9.6的R-loop结果系统性比较(基于K562 cells)
   - 比较DRIP-seq和R-ChIP的peak大小
   - 看两者的peak的overlap
   - 在JUN座位上的R-ChIP, DRIP-seq, DNase, H3K4me1, H3K4me2,H3K4me3, H3K27ac
4. 人类基因组上其他R-loop热点分析
5. R-loop诱导转录
6. R-loop的形成需要游离RNA 末端(free RNA ends)
