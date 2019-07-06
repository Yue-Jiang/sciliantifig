#!/bin/bash 

sample=$1
ref=$2
fq1=$3
fq2=$4
bwa_output_folder=$5

BWA=/net/shendure/vol1/home/yy224/bin/bwa
SAMTOOLS=/net/shendure/vol1/home/yy224/bin/samtools

$BWA mem -M -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" -t 16 $ref $fq1 $fq2 | $SAMTOOLS view -uSb - | $SAMTOOLS sort - > $bwa_output_folder/$sample.srt.bam
$SAMTOOLS index $bwa_output_folder/$sample.srt.bam
