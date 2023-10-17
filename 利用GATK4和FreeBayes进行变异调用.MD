###  计算生物学实验五: exercise instructions for BioHPC Cloud
computers
---


### 1. 实验目的

使用 GATK4 调用短变异（如SNV和InDel）的目的是对基因组数据进行变异分析，以了解个体或种群之间的遗传差异和变异特征。以下是使用 GATK4 调用短变异的一些目的和意义：

1. **突变鉴定**：通过对基因组数据进行短变异调用，可以识别个体或种群中的突变位点，如单核苷酸变异（SNV）和小插入/缺失（Indel）。这有助于了解个体之间的遗传差异，以及与疾病或表型特征相关的变异。
2. **种群遗传学**：通过对不同个体或种群的基因组数据进行短变异调用，可以研究种群之间的遗传多样性、种群结构和进化历史。这有助于了解种群之间的遗传关系、迁移模式和种群适应性。
3. **疾病关联研究**：通过将短变异调用与临床数据或表型特征相关联，可以识别与疾病发生和发展相关的遗传变异。这有助于研究疾病的遗传基础、风险评估和药物治疗的个体化选择。
4. **基因组注释**：通过将短变异调用结果与基因组注释数据库进行比对，可以了解变异位点的功能影响和可能的生物学意义。这有助于解释变异位点对基因功能、调控元件和表达水平的影响。 总之，使用 GATK4 调用短变异可以提供关于个体或种群基因组的详细信息，从而促进遗传学、种群遗传学、疾病研究和基因功能注释等领域的研究。
5. **个体识别和亲子关系分析**：通过分析个体之间的遗传差异，可以用于确定个体的身份，例如在法医学和人口遗传学中。此外，遗传数据还可以用于确定亲子关系，如父子、母子和兄弟姐妹之间的亲缘关系。
6. **变异发生率研究**：通过分析SNV和Indel的分布和频率，可以了解不同类型变异在人群中的发生率。这对于研究遗传多样性和种群历史非常有用。
7. **表观基因组学研究**：通过将短变异调用数据与表观基因组数据（如甲基化）结合使用，可以更全面地理解基因组的调控机制和表达调控。
8. **药物个体化治疗**：通过短变异调用，可以识别与药物代谢、反应和不良事件相关的遗传变异，从而帮助医生为患者选择最合适的治疗方案，实现个体化医疗。
9. **基因功能研究**：研究SNV和Indel的分布和功能影响有助于理解基因功能和调控机制，进一步推进基因组学和生物学领域的研究。
10. **突变研究**：通过调用短变异，可以发现新的突变位点，这对于疾病致病突变的鉴定以及研究基因组变异的新机制具有重要意义。

### 2. 实验准备

#### 2.1 实验平台

Linux JSvr01 3.10.0-1160.80.1.el7.x86_64 #1 SMP Tue Nov 8 15:48:59 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

#### 2.2 数据简述

| 数据                          | 来源        |
| ----------------------------- | ----------- |
| genome.fa                     | Linux服务器 |
| SRR1663608_thinned_1.fastq.gz | Linux服务器 |
| SRR1663608_thinned_2.fastq.gz | Linux服务器 |
| SRR1663609_thinned_1.fastq.gz | Linux服务器 |
| SRR1663609_thinned_2.fastq.gz | Linux服务器 |
| SRR1663610_thinned_1.fastq.gz | Linux服务器 |
| SRR1663610_thinned_2.fastq.gz | Linux服务器 |
| SRR1663611_thinned_1.fastq.gz | Linux服务器 |
| SRR1663611_thinned_2.fastq.gz | Linux服务器 |



#### 2.3 软件配置

| 软件     | 系统  | 版本         |
| -------- | ----- | ------------ |
| bwa      | Linux | 0.7.17-r1188 |
| GATK     | Linux | 4.1.9.0      |
| Samtools | Linux | 1.7          |
| vcftools | Linux | 0.1.15       |



### 3. 实验内容

我们将使用黑腹果蝇 WGS 双端 Illumina 数据和 NCBI 登记号 SRR1663608、SRR1663609、SRR1663610、SRR1663611，分别对应于样品 ZW155、ZW177、ZW184 和 ZW185。为了加快计算速度，最初每个样本约 1000 万个读取对的数据已被下采样 50%。

该练习包括从 WGS Illumina 读取到四个样本变异调用的几个步骤，如 GATK 最佳实践中所指定。要调用所有四个样本中的变体，需要对四个 FASTQ 文件对中的每一个执行从读取比对到单倍型调用的所有步骤。

在具有大内存的多 CPU 计算机上，可以并行启动此类运行（例如，通过在后台手动提交给定步骤，每次针对不同的样本）。

为了练习的目的，请针对您选择的至少一个样本运行所有这些步骤。

最终联合变体调用步骤所需的其他样本的所有中间结果都已预先计算，并且可以从研讨会目录中复制（请参阅下面的详细说明）。

**注意：这个部分不可以为空*

3.1准备工作目录

创建目录并将 /shared_data/annotation2018/ 中的数据文件复制到 /home/S2021281017/lab5中，。您还将创客软件目录复制到/workdir/USER。制造商软件目录包括大型序列重复数据库。

```bash
mkdir lab5mkdir lab5mkdir lab5mkdir lab5mkdir lab5mkdir lab5
cd lab4
cp /work/cb-data/Variants_workshop_2020/*.fastq.gz .
cp -r /work/cb-data/Variants_workshop_2020/scripts .
cp -r /work/cb-data/Variants_workshop_2020/genome .
ls -al
```

![image-20231016112302772](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016112302772.png)

> 以上是准备好的文件

检查文件大小 – 每个文件大小应约为 410-430 MB。除了读取的文件之外，子目录 script 中还提供了几个 shell 脚本 *.sh，对应于管道的各个阶段。子目录基因组包含一个带有黑腹果蝇参考基因组的 FASTA 文件基因组.fa。

#### 3.2准备参考基因组

子目录genome中的文件genome.fa是我们进行对比的参考。在启动管道之前，需要对基因组进行BWA比对索引。同时，需要对所有染色体长度的信息进行汇总。生成了在子目录genome下的一系列文件。文件genome.fa.fai和genome.fa.dict是简单的文本文件，总结了染色体大小和原始FASTA文件中的起始字节位置。其他文件构成BWA索引。

所有这些均通过运行shell脚本prepare_genome.sh来完成

```bash
cd scripts
chmod u+x prepare_genome.sh
vim prepare_genome.sh
```

> 修改脚本权限和内容，使其能够被执行

![image-20231016145550771](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016145550771.png)

> 这个脚本是用于准备一个参考基因组，以便进行生物信息学分析，特别是在变异调用和对比中使用，下面逐步解释脚本中每个部分
>
> 'TMP' = ...这一行定义了名为TMP的变量，其中包含了特定的目标路径。这是用于在脚本执行期间存储临时文件的目录
>
> GATKDIR =  . . .这一行使GATKDIR选择好了目标路径，GATK是基因组工具分析包
>
> export PATH这一行将GATKDIR添加到系统的PATH环境变量中
>
> 'cd genome'将当前的工作目录更改为名为‘genome'的子目录。

```bash
cd /home/S20213281017/lab5
nohup ./scripts/prepare_genome.sh >& prepare_genome.log &
top -u S20213281017
cd genome
ls -al
```

该程序将在后台运行，将所有屏幕输出保存到命名为 bwa_aln_SRR1663609.log 的日志文件中。大约运行时间：10 分钟。可以用top -u来监视进程，完成后，列出子目录genome.fa.fai 和 genome.fa.dict是简单的文本文件，总结了染色体大小和原始FASTA文件中的起始位置。其他文件构成BWA索引。

结果将是文件 SRR1663609.bam ，以 BAM 格式描述比对。该文件是一个二进制文件。

![image-20231016152356958](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016152356958.png)

> 如图所示，得到genome.fa.fai等文件

```bash
vim prepare_genome.log
```

![image-20231016153512772](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016153512772.png)

> 查看prepare_genome.log

```bash
vim genome.fa.fai
```

![image-20231016153804590](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016153804590.png)

> 第一列为染色体名称，第二列为染色体长度，第三列为偏移量，第四列为每行的碱基数，第五列为每行的字符数，包括换行符

```bash
vim genome.dict
```

![image-20231016212531572](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016212531572.png)

> ##### @SQ用于分割每一条染色体信息，第一列是染色体信息，第二列是染色体名字，第三列是M5码，第四列是来源

#### 3.3排序并标记重复项

在此步骤中，之前获得的 BAM 文件中发现的任何重复片段现在将被标记为在下游分析中忽略，并且比对将根据基因组坐标进行排序。排序后的 BAM 文件也将被索引。所有这些操作现在均由名为 MarkDuplicatesSpark 的单个 GATK 工具执行。

##### 3.3.1 将样本与参考基因组比对

首先修改写好的脚本权限并启动脚本

```bash
cd scripts
chmod u+x bwa_aln.sh
cd ..
nohup ./scripts/bwa_aln.sh SRR1663608 ZW155 >& bwa_aln_SRR1663608.log &
nohup ./scripts/bwa_aln.sh SRR1663609 ZW177 >& bwa_aln_SRR1663609.log &
nohup ./scripts/bwa_aln.sh SRR1663610 ZW184 >& bwa_aln_SRR1663610.log &
nohup ./scripts/bwa_aln.sh SRR1663611 ZW185 >& bwa_aln_SRR1663611.log &
ls -al
```

![image-20231016220750394](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016220750394.png)

> 得到如上图所示四个bam文件

bam文件是一个二进制文件，利用samtools命令可以检查它的内容，它读取这个二进制文件并将其转换为人类可读的文本格式，打印到屏幕上：

```bash
samtools view -h SRR1663608.bam | less
```

![image-20231016221332018](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231016221332018.png)

> 以SRR1663609.bam文件为例
>
> 使用Samtools查看 BAM 文件后，每行都代表比对的一个测序读取。SAM（Sequence Alignment/Map）格式和BAM格式都包含类似的信息。以下是 SAM 文件或 BAM 文件中各列的含义：
>
> 1. **QNAME（列1）**：Query名，它是每个测序读取的唯一标识符。
> 2. **FLAG（列2）**：标志位，指示比对的属性，例如是否匹配、是否为第一片段、是否为最后一片段等。此列的值通常以十进制数表示，但也可以使用十六进制来表示。
> 3. **RNAME（列3）**：Reference名，这是比对到的参考基因组的名称。如果无法比对到参考基因组，则该字段通常为"*"。
> 4. **POS（列4）**：比对的起始位置，表示测序读取在参考基因组上的起始位置。它是以1为起点的，即第一个碱基的位置为1。
> 5. **MAPQ（列5）**：比对质量，它表示比对的可信度。通常以Phred质量分数的形式给出，数值越高表示比对质量越好。
> 6. **CIGAR（列6）**：比对信息字符串，它描述了测序读取与参考基因组的比对方式，包括插入、删除、匹配等信息。
> 7. **RNEXT（列7）**：下一个比对的Reference名，用于标识成对的测序读取。如果没有成对的测序读取或无法比对到参考基因组，则该字段通常为"*"。
> 8. **PNEXT（列8）**：下一个比对的起始位置，用于标识成对的测序读取的起始位置。
> 9. **TLEN（列9）**：测序读取的长度，通常是负值，表示测序读取成对的片段之间的距离。
> 10. **SEQ（列10）**：测序读取的碱基序列。
> 11. **QUAL（列11）**：测序质量，表示测序读取每个碱基的质量得分。通常以Phred质量分数的形式给出。
>
> 可以看出，这一步完成了两组信息的比对，为下一步对重复的reads进行排序和标记做准备

##### 3.3.2对重复reads进行排序和标记

在这一步中，之前获得的BAM文件中发现的重复片段被标记，并在下游分析中被忽略，再根据基因组坐标对比进行排序。结果排序的BAM文件将被索引。所有的这些操作都由一个名为MarkDuplicateSpark的GATK工具执行.依照sort_dedup_index.sh中的使用，最终会产生sorted.dedup.bam文件，以及索引文件sortrd.dedup.bai文件

首先修改脚本权限并执行脚本

```bash
cd scripts
chmod u+x sort_dedup_index.sh
vim sort_dedup_index.sh #打开脚本修改
```

![image-20231017001926181](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017001926181.png)

> 将TMP所指路径换成工作路径即可

```bash
nohup ./scripts/sort_dedup_index.sh SRR1663608 >& sort_dedup_index_SRR1663608.log &
nohup ./scripts/sort_dedup_index.sh SRR1663609 >& sort_dedup_index_SRR1663609.log &
nohup ./scripts/sort_dedup_index.sh SRR1663610 >& sort_dedup_index_SRR1663610.log &
nohup ./scripts/sort_dedup_index.sh SRR1663611 >& sort_dedup_index_SRR1663611.log &
```

![image-20231017005655305](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017005655305.png)

运行 samtools 命令检查获得的文件的对齐统计摘要,以其中一个为实例

```bash
samtools flagstat SRR1663609.sorted.dedup.bam
```

![image-20231017010151006](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017010151006.png)

> 重复reads的数量为750524,而通过QC的reads数量为10209495
>
> 可以看出比对到参考基因组上的reads数量为9698015，占比94.99%
>
> 正确地匹配到参考序列的reads数量为8749418，占比86.3%

#### 3.4用IGV可视化对齐结果

利用Xftp将SRR1663609.sorted.dedup.bam、SRR1663609.sorted.dedup.bam.bai以及genome.fa传输至电脑桌面，然后用IGV可视化

![image-20231017013050334](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017013050334.png)

> 在IGV图中，灰色表示insert size正常，是常规reads。 红色表示insert size大于预期，表示可能存在缺失突变。 蓝色表示insert size小于预期，表示可能存在插入突变。而其他颜色的reads表示paired end reads中成对两条reads比对到了不同的染色体。比如在Chr2R:826964-827007的配对的reads对比到了chr3RHet:1819995上了。除此之外，紫色的部分说明了存在插入碱基，数字表示插入碱基的个数。当覆盖轨道为彩色时，表示核苷酸与参照碱基比对的差异查过20%,可以根据这个，判断这些点位的SNP。

#### 3.5在单个样本上运行GATK HaplotypeCaller

GATK（Genome Analysis Toolkit）的HaplotypeCaller是一个用于变异调用的工具，特别是用于单核苷酸多态性（SNP）和插入/删除多态性（Indel）的检测。它的主要功能包括：

1. **局部重组（Local Reassembly）**：HaplotypeCaller采用了一种先进的方法，对重组区域进行局部重组。这意味着它能够捕获到复杂的变异，如Indel或者SNP，尤其在重组区域，提高了变异检测的准确性。
2. **不考虑PCR复制（PCR Indel Model）**：该工具在变异检测时，考虑到了PCR扩增引入的伪变异问题，从而提高了变异检测的准确性。
3. **多样性调用（Allele-Specific Annotation）**：HaplotypeCaller不仅返回了每个位点的变异信息，还提供了多样性信息，如每个变异的可能等位基因和每个样本的各个等位基因的深度信息。这对于进一步的种群遗传学和个体遗传学分析非常有用。
4. **标定质量分数（Recalibration）**：HaplotypeCaller还支持GATK中的质量分数重新校准，以提高变异检测的质量。

在此步骤中，我们将在 BAM 文件上运行 HaplotypeCaller，以生成感兴趣基因组区域中每个基因座的每个样本的基因型可能性信息。通常，该区域是整个基因组。为了节省时间，我们将专注于一条染色体，即chr2R。

首先修改脚本内容中TMP位置以及修改脚本权限

```bash
cd scripts
vim u+x hc.sh
chmod u+x hc.sh
```

![image-20231017014404180](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017014404180.png)

> 已修改的shell文件

运行脚本

```bash
nohup ./scripts/hc.sh SRR1663608 >& hc_SRR1663608.log &
nohup ./scripts/hc.sh SRR1663609 >& hc_SRR1663609.log &
nohup ./scripts/hc.sh SRR1663610 >& hc_SRR1663610.log &
nohup ./scripts/hc.sh SRR1663611 >& hc_SRR1663611.log &
```

```bash
ls
```

得到运行得到的文件，应是四个g.vcf文件，除此以外还有日志文件

![image-20231017141223886](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017141223886.png)

查看.g.vcf文件

```bash
vim SRR1663609.g.vcf
```

![image-20231017140438750](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017140438750.png)

> 如图所示为前面的注释部分

![image-20231017140719205](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017140719205.png)

> 主体部分包含十列数据，主体部分每一行代表一个变异位点信息
>
> VCF（Variant Call Format）文件是一种用于描述遗传变异数据的标准文件格式。每一列都有特定的含义，下面是每一列的解释：
>
> 1. **#CHROM**：这是染色体的名称或标识符，通常是染色体的编号或名称。可以是字符串或数字。
> 2. **POS**：这是变异的位置，通常是相对于染色体的位置。
> 3. **ID**：这是变异的唯一标识符，通常是一个rs号或其他唯一的标识符。如果没有可用的唯一标识符，通常是"."。
> 4. **REF**：这是参考基因组中的碱基。在该位置，参考基因组的碱基。
> 5. **ALT**：这是代表变异的碱基或碱基序列。可能有多个不同的变异类型，因此这一列可以包含多个不同的值，用逗号分隔。
> 6. **QUAL**：这是一个质量分数，表示这个变异的质量。通常是一个浮点数，用来表示可信度或置信度。
> 7. **FILTER**：这是一个筛选器字段，指示该变异是否通过了一些质量控制筛选。可以是PASS，或者是一个代表筛选不过的代码或标识符。
> 8. **INFO**：这是一个包含关于变异的额外信息的列，通常是以分号分隔的键值对。
> 9. **FORMAT**：这一列定义了后续列中样本的数据格式。它通常包含了一个或多个标签，用于描述后续列中的数据。

#### 3.6使用GenotypeGVCFs进行联合变异的调用

##### 3.6.1首先将四个vcf文件合并

在中间样本级文件 *.g.vcf 结果可用于联合调用所有四个样本上的变体之前，必须使用 GATK4 的工具 JointGVCFs 将它们合并为单个多样本 g.vcf 文件。调用该命令的脚本combineGVCFs.sh

```bash
cd scripts
chmod u+x combineGVCFs.sh
chmod u+x genotypeGVCFs.sh
vim combineGVCFs.sh
```

![image-20231017142710860](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017142710860.png)

> 修改tmp到工作目录

```bash
cd ..
nohup ./scripts/combineGVCFs.sh >& combineGVCFs.log &
ls
```

![image-20231017143155327](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017143155327.png)

> 如图所示生成了all.g.vcf，即合并成功

```bash
cd scripts
vim genotypeGVCFs.sh
```

![image-20231017144319312](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017144319312.png)

> 修改tmp到工作路径

```bash
cd ..
nohup ./scripts/genotypeGVCFs.sh >& genotypeGVCFs.log &
ls
```

![image-20231017144744250](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017144744250.png)

> 如图所示生成了all.vcf文件

#### 3.7使用FreeBayes进行联合变异调用

FreeeBayes 是 Erik Garrison 等人开发的变体调用程序，https://github.com/ekg/freebayes。它是 GATK 的替代品之一。与 GATK 的 HaplotypeCaller 类似，FreeBayes 使用基于单倍型的方法进行变异检测，尽管实现方式不同。FreeeBayes 的输入包含所有涉及样本的比对 BAM 文件。我们将使用本练习之前获得的经过处理的 BAM 文件

首先修改脚本内容和权限

```bash
cd scripts
chmod u+x fb.sh
vim fb.sh
```

![image-20231017145415908](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017145415908.png)

> 将FB路径修改到目的路径

```bash
cd ..
nohup ./scripts/fb.sh &> fb.log
ls
```

![image-20231017150705149](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017150705149.png)

> 生成了fb.vcf文件

```bash
vim fb.vcf
```

![image-20231017150903522](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017150903522.png)

> 如图所示，与前面的vcf文件类似

#### 3.8过滤变异体

前面步骤中获得的 VCF 文件是原始结果，可能包含大量误报，具体取决于调用中使用的严格选项。由于调用步骤非常耗时，因此通常建议设置这些选项以发出一组包含的变体，然后通过各种参数（例如 VCF 文件的 INFO 字段中记录的参数）过滤该组。

##### 3.8.1利用grep工具过滤

在GATK中选项 --stand-call-conf 控制要输出的变体质量（VCF 的 QUAL 字段）的下限。该选项应设置为某个较低的值（例如 5，而不是默认的 30）。可以使用许多不同的工具来完成对原始变体集的过滤。在很多情况下，人们可以简单地使用标准的 Linux 文本解析工具，例如 grep、awk 或 sed。例如，要从原始 all.vcf 中提取 QUAL 大于（例如 60）的变体子集，我们可以使用以下命令：

```bash
grep '#' all.vcf > all.qua160.vcf
vim all.qua160.vcf
```

![image-20231017154308460](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017154308460.png)

> 将 VCF 标题行（包含“#”）提取到新的 VCF 文件中

```bash
grep -v "#" all.vcf | awk '{if($6>60) print}' >> all.qual60.vcf
vim all.qual60.vcf
```

![image-20231017155539576](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017155539576.png)

> 处理非标题行，仅当第六列（即 QUAL 所在的位置）高于 60 时才将它们附加到新文件中

##### 3.8.2使用VariantFiltration过滤



GATK 提供了一个名为 VariantFiltration 的工具，它允许更复杂的过滤模式。使用此工具的示例脚本称为 filter_vcf.sh 。使用此脚本时，应用于 SNP 的过滤标准与应用于 indel 的过滤标准不同。为了实现此目的，首先将 SNP 和 indel 提取到单独的文件中，对这些文件进行过滤，然后将 SNP 和 indel 过滤后的文件合并回单个过滤文件中。

首先修改脚本权限及脚本权限

```bash
cd scripts
chmod u+x filter_vcf.sh
vim filter_vcf.sh
```

![image-20231017160509094](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017160509094.png)

> 修改tep路径到工作路径

```bash
cd ..
nohup ./scripts/filter_vcf.sh all >& filter_vcf.log &
ls
```

![image-20231017161101732](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017161101732.png)

> 如图所示，生成all.filtered.vcf文件

（请注意，我们提供 VCF 文件名的前缀，即不带 .vcf 扩展名作为参数）。过滤后的VCF文件将被称为all.filtered.vcf（也会创建相应的索引文件*.vcf.idx）。其他中间文件（带有单独的 SNP、插入缺失、过滤的和未过滤的以及它们的索引）也将被生成 - 这些文件可能会被删除。

检查过滤后的 VCF 文件。请注意 FILTER 字段中的变化。而不是点“.” （没有过滤信息），该字段现在将包含标志 PASS（通过过滤器的变体）和 my_snp_filter 或 my_indel_filter（这两个字符串都是在过滤命令中定义的）——标记未通过相应过滤器的变体。请注意，没有从文件中删除任何变体。过滤失败的那些只是被标记为这样。

```bash
vim all.filtered.vcf
```

![image-20231017161340723](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017161340723.png)

#### 3.9基本统计数据和变体集比较

##### 3.9.1使用Linux文本解析工具处理

给定一个 VCF 文件，其最简单的属性可以通过运行标准 Linux 文本解析工具来获得

例如，要获取文件中的变量数量，请运行以下命令：

```bash
grep -v "#" all.vcf | wc -l
```

![image-20231017162116342](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017162116342.png)

> grep 过滤掉标题行并将其输出传输到 wc -l 中，后者计算剩余行并将结果显示在屏幕上

要提取位于染色体 chr2R 上位置 10000 和 20000 之间的位点并将其保存在文件中，只需运行

```bash
grep -v "#" all.vcf | awk '{if($1=="chr2R" && $2 >=10000 && $2 <=20000) print}' > extracted_records
vim extracted_records
```

![image-20231017162627853](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017162627853.png)

> 请注意，在这种情况下，染色体条件 $1=="chr2R" 并不是真正需要的，因为我们的 VCF 文件仅包含 chr2R 的数据，但是，对于更通用的输入来说，这是必要的。

要快速查找有多少变量通过了过滤需要以下命令

```bash
awk '{if($7=="PASS") print}' all.filtered.vcf | wc -l
```

![image-20231017163028942](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017163028942.png)

> 得到有417825通过了过滤

对 VCF 文件更复杂的分析和操作可以使用专门的软件工具来完成，例如 GATK 包中包含的工具和 vcftools 包（独立于 GATK）中的工具。

##### 3.9.2使用 GATK4 的一致性和基因型一致性工具

GATK4 提供了有趣的函数 Concordance 和 GenotypeConcordance，用于总结变体集的各种统计数据，并将其与从相同数据但使用不同方法获得的另一个变体集进行比较。基于此函数的脚本 var_compar.sh 将比较任意两个 VCF 文件，例如 all.vcf 和 fb.vcf（分别使用 GATK4 和 FreeeBayes 获得）：

首先修改脚本的内容和权限

```bash
cd scripts
chmod u+x var_compar.sh
vim var_compar.sh
```

![image-20231017163759518](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017163759518.png)

> 修改tmp到工作路径

```bash
cd ..
./scripts/var_compar.sh fb all ZW155 >& var_compar.log &
ls
```

![image-20231017164555617](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017164555617.png)

> 可以看到得到了fb.all.comp.site_summary这些文件

```bash
vim fb.all.comp.site_summary
```

![image-20231017164936881](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017164936881.png)

> 以all.vcf为原假设，Precition = 0.792

##### 3.9.3使用vcftools工具

vcftools(A.Auton，A.Amrcketta，https://vcftools.github.io/index.html))是分析和操作vcf文件的流行工具包。

获取VCF基本统计数据(样本数和变异站点数)：

```bash
vcftools --vcf all.vcf
```

![image-20231017165545700](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017165545700.png)

> 共有四个样本，变异点数为423493

提取变异子集(染色体chr2R，位置1M和2M之间)并将其写入新的VCF文件

```bash
vcftools --vcf all.vcf --chr chr2R --from-bp 1000000 --to-bp 2000000 --recode --recode-INFO-all -c > subset.vcf
vim subset.vcf
```

![image-20231017170224480](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017170224480.png)

获取所有变体的等位基因频率并将其写入文件

```bash
vcftools --vcf all.vcf --freq -c > all.freqs
vim all.freqs
```

![image-20231017170446143](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017170446143.png)

比较两个 VCF 文件（将在文件 fb.all.compare.* 中打印出各种比较信息）

```bash
vcftools --vcf fb.vcf --diff all.vcf --out fb.all.compare
```

![image-20231017170705044](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20231017170705044.png)

发现fb.vcf和all.vcf，都至少含有399691个短变异，但从以上实验可知，使用GATK4和FreeBayes进行联合变异的调用，GATK4能调用更多的短变异体。

### 4. 实验总结

#### 4.1 实验结论

GATK4（Genome Analysis Toolkit）和FreeBayes都是用于变异调用的工具，但它们在算法和方法上有一些重要的区别。以下是它们之间的主要区别：

1. 数据处理方法：
   - GATK4：GATK是一套广泛使用的工具，旨在减少测序数据的误差，包括基础质量校正、局部再比对等步骤，以提高变异检测的准确性。它还使用多个不同的过滤器来过滤潜在的假阳性突变。GATK4使用贝叶斯方法来估计突变的概率。
   - FreeBayes：FreeBayes是一个基于贝叶斯统计的工具，它直接从原始测序数据中估计变异的概率，而不需要复杂的数据处理步骤。它不执行像GATK那样的多样的数据处理，而是更专注于从原始数据中获取变异。
2. 引入数据的方式：
   - GATK4：GATK通常需要一个已知的参考基因组，然后将测序数据比对到该参考基因组上。它的变异检测过程通常包括从比对的测序数据中识别SNPs（单核苷酸多态性）和Indels（插入和缺失）。
   - FreeBayes：FreeBayes也通常需要一个参考基因组，但它可以使用比对或对齐的方式，也可以使用未比对的测序数据。
3. 调用策略：
   - GATK4：GATK在调用变异时使用多种统计模型，包括单样本和多样本联合分析。它使用一种叫做HaplotypeCaller的工具来调用变异，以获得更好的变异分辨率。
   - FreeBayes：FreeBayes使用贝叶斯模型来进行变异调用，通常是在每个样本独立进行。虽然也可以进行多样本分析，但它的主要焦点是单样本变异调用。

#### 4.2 实验收获

在本次实验中主要使用到了GATK4这个强大的生物信息学工具，下面是用其调用短变异体的过程总结：

1. 数据预处理：
   - 数据质量控制：首先，您需要执行数据质量控制步骤，包括去除低质量的读数、去除接头序列、修复碱基质量分数等。这可以使用工具如FastQC和Trimmomatic来完成。
2. 比对到参考基因组：
   - 用BWA、Bowtie2等比对工具将测序数据比对到参考基因组。这将创建一个BAM文件，其中包含了每个读数的比对信息。
3. 数据处理和改进：
   - GATK推荐使用工具如Picard或Samtools来处理BAM文件，执行重复标记、本地重比对、碱基质量校正等步骤。这有助于提高变异检测的准确性。
4. 比对质量再校正（BQSR）：
   - 使用GATK的BaseRecalibrator工具来估计碱基质量分数的校正表。这有助于减少碱基质量对变异调用的影响。
5. 变异调用：
   - 使用GATK的HaplotypeCaller工具来进行变异调用。HaplotypeCaller采用了一种全局模型，它会考虑整个区域的信息，而不是单个读数。这有助于提高SNP和Indel的检测准确性。
6. 变异筛选和过滤：
   - GATK提供了多个过滤和筛选工具，您可以根据不同的标准来选择筛选标准，以过滤掉可能的假阳性变异。
7. 多样本联合分析（可选）：
   - 如果有多个样本，可以使用GATK的GenotypeGVCFs工具来进行多样本联合分析，从而更好地识别变异。
8. 结果报告：
   - 最终，可以使用GATK的VariantsToTable或VariantsToVCF工具将变异调用结果导出到相应的文件格式，以供后续分析和注释使用。



#### 4.3 其它心得

这个不是必须项目
