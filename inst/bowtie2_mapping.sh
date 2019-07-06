#!/bin/bash 

sample=$1
ref=$2
fq1=$3
fq2=$4
bt2_output_folder=$5

BOWTIE2=/net/shendure/vol1/home/yy224/src/miniconda2/bin/bowtie2
SAMTOOLS=/net/shendure/vol1/home/yy224/bin/samtools
$BOWTIE2 -p 8 -X 2000 -3 1 -x $ref -1 $fq1 -2 $fq2 | $SAMTOOLS view -uSb - | $SAMTOOLS sort - > $bt2_output_folder/$sample.srt.bam
$SAMTOOLS index $bt2_output_folder/$sample.srt.bam

