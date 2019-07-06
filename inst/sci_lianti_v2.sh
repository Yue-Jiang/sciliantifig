#!/bin/bash
#sci-ligation-lianti version 2 pipeline

fastq_folder=$1 # the folder for fastq files
all_output_folder=$2 # the output folder
sample_ID=$3 # the sample ID for each PCR samples after demultiplex
core=$4 # core number for computation
cutoff=$5  # the number of unique reads cutoff for splitting single cells
barcodes_tn5=$6 # the round1 barcode list (on tn5) for splitting single cells
barcodes_ligation=$7 # the round2 barcode list (by ligation) for splitting single cells
barcodes_sss=$8 # the SSS barcode list for splitting single cells
script_folder=$9 # the script folder for called python scripts

#define the mismatch rate (edit distance) of UMIs for removing duplicates:
mismatch_tn5=1
mismatch_ligation=1
mismatch_sss=0
mismatch_RT=3
RT_primer="GGGATGCAGCTCGCTCCTG"

#define the bin of python
python_path="/net/shendure/vol1/home/yy224/src/miniconda2/bin/"

#define the path to genome reference
# ref_bwa="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m.fa"
# ref_bt2="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m"
ref_bwa="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/mm10.fa"
ref_bt2="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/mm10_natural_order"
# ref_bwa="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/Spretus.pseudo.fa"
# ref_bt2="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/Spretus.pseudo"
# ref_bwa="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/dm3.fa"
# ref_bt2="/net/shendure/vol8/projects/fly_recomb/lh3/bwadb/dm3"

#define the location of script:
script_path=$script_folder

SAMTOOLS=/net/shendure/vol1/home/yy224/bin/samtools
BEDTOOLS=/net/shendure/vol1/home/yy224/bin/bedtools
CUTADAPT=/net/shendure/vol1/home/yy224/src/miniconda2/bin/cutadapt

############ split by SSS barcode
#get SSS barcode with specified mismatch, evaluate percentage of junk sequence from the truseq prep
###############################################################
#######################step 1, up to step2
###############################################################
# input_folder=$fastq_folder
# output_folder=$all_output_folder/split_SSS
# script=$script_folder/order_read_pairs.py
# echo "changing the name of fastq files..."
# 
# for sample in $(cat $sample_ID); do echo changing name $sample; ln -s $input_folder/$sample\_*R1*gz $input_folder/$sample.R1.fastq.gz; ln -s $input_folder/$sample\_*R2*gz $input_folder/$sample.R2.fastq.gz; done
# 
# echo "Attaching sss barcode"
# mkdir -p $output_folder
# $python_path/python $script $input_folder $sample_ID $output_folder $barcodes_sss $core $mismatch_sss $mismatch_RT $RT_primer
# echo "SSS barcode attached."

###############################################################
#######################step 2, up to step 3
###############################################################
# input_folder=$fastq_folder
# output_folder=$all_output_folder/split_SSS
# script=$script_folder/order_read_pairs.py
# 
# # mkdir -p $output_folder/bkp
# # cp $output_folder/*.ordered.* $output_folder/bkp
# # After generating and backing up the "ordered" file, change sample_ID file to only include lib121, concatenate nextseq171010 and nextseq171017 yi117 sequences together, run the split by SSS code below on the merged files in the 20171010 folder
# 
# echo "split by sss"
# cat $barcodes_sss | while read line
# do
# cat $sample_ID | while read lib
# do
# echo "zcat $output_folder/$lib.R1.ordered.fastq.gz | grep ^@$line -A3 | grep -v -- ^--$ | gzip - > $output_folder/$lib.R1.$line.fq.gz" | qsub -cwd -V -l mfree=2G -N splitR1sss
# echo "zcat $output_folder/$lib.R2.ordered.fastq.gz | grep ^@$line -A3 | grep -v -- ^--$ | gzip - > $output_folder/$lib.R2.$line.fq.gz" | qsub -cwd -V -l mfree=2G -N splitR2sss
# done
# done
# echo "split by sss completed"
 
###############################################################
#######################step 3, up to step 4
#######################prepare $all_output_folder/split_sss_sample_sheet.txt file before running
#######################reformat $output_folder/sss_summary.txt afterwards 
#######################remove sss barcode not used afterwards
###############################################################
# input_folder=$fastq_folder
# output_folder=$all_output_folder/split_SSS
# script=$script_folder/order_read_pairs.py
# 
# echo "get sss summaries"
# cat $barcodes_sss | while read line
# do
# cat $sample_ID | while read lib
# do
# echo $lib $line >> $output_folder/sss_barcode_count.txt
# zcat $output_folder/$lib.R1.$line.fq.gz | wc -l  | awk '{print $1/4}' >> $output_folder/sss_barcode_count.txt
# done
# done
# 
# cat $all_output_folder/split_sss_sample_sheet.txt | while read -r r1 r2 r3 r4
# do
# echo $r1 $r2 $r3 $r4 >> $output_folder/sss_summary.txt
# grep "$r1 $r4" -A1 $output_folder/sss_barcode_count.txt >> $output_folder/sss_summary.txt
# done
# echo "sss summaries completed"

###############################################################
#######################step 4, up to step5 
###############################################################
############ split by tn5 and ligation barcodes
####get tn5 and ligation barcodes with specified mismatch, evaluate percentage of junk sequence from tagmentation and ligation
# echo
# echo "Start trimming the read1 file..."
# echo $(date)
# split_SSS_folder=$all_output_folder/split_SSS
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# script=$script_folder/parse_tn5_ligation_barcodes.py
# sample_tn5_ligation_file=$trimmed_fastq_folder/sample_tn5_ligation.txt
# 
# mkdir -p $trimmed_fastq_folder
# 
# echo "Cutadapt R1..."
# cat $split_SSS_folder/sss_summary.txt | while read -r lib NA cell_number sss read_number
# do
# sem -j $core $CUTADAPT -g AGATGTGTATAAGAGACAG -O 19 -e 0.2 --info-file $trimmed_fastq_folder/$lib.R1.$sss.info.tab -o $trimmed_fastq_folder/$lib.R1.$sss.trimmed.fq.gz $split_SSS_folder/$lib.R1.$sss.fq.gz
# #$CUTADAPT -g AGATGTGTATAAGAGACAG -O 19 -e 0.2 --info-file $trimmed_fastq_folder/$lib.R1.$sss.info.tab -o $trimmed_fastq_folder/$lib.R1.$sss.trimmed.fq.gz $split_SSS_folder/$lib.R1.$sss.fq.gz
# done
# sem --semaphoretimeout 1800
# echo "All R1 trimmed file generated."
# 
# wc -l $split_SSS_folder/*.info.tab

###############################################################
#####################wc -l *.info.tab to confirm###############
#######################step 5, up to step6 
###############################################################
# split_SSS_folder=$all_output_folder/split_SSS
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# script=$script_folder/parse_tn5_ligation_barcodes.py
# sample_tn5_ligation_file=$trimmed_fastq_folder/sample_tn5_ligation.txt
# 
# echo "reformatting info file"
# cat $split_SSS_folder/sss_summary.txt | while read -r lib NA cell_number sss read_number
# do
# awk -F $'\t' '{if($2>=0) {print substr($5,length($5)-20,7), substr($5,length($5)-7,8);} else print "NA","NA"}' $trimmed_fastq_folder/$lib.R1.$sss.info.tab  | gzip - > $trimmed_fastq_folder/$lib\_$sss.R1.bc1.bc2.txt.gz
# done
# echo "info file reformatted"
# mkdir -p $trimmed_fastq_folder/bkp
# mv $trimmed_fastq_folder/*.info.tab $trimmed_fastq_folder/bkp
# #gzip $trimmed_fastq_folder/bkp/*.info.tab
# 
# cat $split_SSS_folder/sss_summary.txt | while read -r lib NA cell_number sss read_number
# do
# ln -s $trimmed_fastq_folder/$lib.R1.$sss.trimmed.fq.gz $trimmed_fastq_folder/$lib\_$sss.R1.fq.gz
# ln -s $split_SSS_folder/$lib.R2.$sss.fq.gz $trimmed_fastq_folder/$lib\_$sss.R2.fq.gz
# echo "${lib}_${sss}" >> $trimmed_fastq_folder/sample_tn5_ligation.txt
# done
# 
# echo "Attaching tn5 and ligation barcodes"
# $python_path/python $script $trimmed_fastq_folder $sample_tn5_ligation_file $barcodes_tn5 $barcodes_ligation $core $mismatch_tn5 $mismatch_ligation
# echo "Tn5 and ligation barcodes attached"

###############################################################
#######################step 5, up to step6 
#######################wc -l $trimmed_fastq_folder/sample_tn5_ligation.txt
###############################################################
# split_SSS_folder=$all_output_folder/split_SSS
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# script=$script_folder/parse_tn5_ligation_barcodes.py
# sample_tn5_ligation_file=$trimmed_fastq_folder/sample_tn5_ligation.txt
# 
# echo "trimming R2"
# cat $sample_tn5_ligation_file | while read line
# do
# sem -j $core $CUTADAPT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -O 13 -e 0.2 -o $trimmed_fastq_folder/$line.R2.trimmed.attached.fastq.gz $trimmed_fastq_folder/$line.R2.attached.fastq.gz
# # $CUTADAPT -g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -O 13 -e 0.2 -o $trimmed_fastq_folder/$line.R2.trimmed.attached.fastq.gz $trimmed_fastq_folder/$line.R2.attached.fastq.gz
# done
# sem --semaphoretimeout 1800
# echo "All trimmed file generated."
# 
###############################################################
#######################step 6, up to here, 
#######################
###############################################################

#mv $trimmed_fastq_folder/*.R2.attached.* bkp
#mv $trimmed_fastq_folder/*.trimmed.fq.gz bkp
#find . -name "*" -size -200c -delete #delete the ln -s files

############ align by bwa mem
##define the output folder
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# processed_fastq_folder=$trimmed_fastq_folder
# bwa_output_folder=$all_output_folder/spret_bwa_alignment
# processed_fastq_list=$trimmed_fastq_folder/sample_tn5_ligation.txt
# script=$script_folder/bwa_mem.sh
# 
# echo "Start alignment using bwa mem..."
# echo input folder: $processed_fastq_folder
# echo processed_fastq_list: $processed_fastq_list
# echo output_folder: $bwa_output_folder
# mkdir -p $bwa_output_folder
# #start the alignment
# for sample in $(cat $processed_fastq_list); do echo Aligning $sample; echo "$script $sample $ref_bwa $processed_fastq_folder/$sample.R1.trimmed.attached.fastq.gz $processed_fastq_folder/$sample.R2.trimmed.attached.fastq.gz $bwa_output_folder" | qsub -cwd -V -l mfree=20G -N mem; done
# echo "All mem alignment done."

########### align by bowtie2
#define the output folder
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# processed_fastq_folder=$trimmed_fastq_folder
# bt2_output_folder=$all_output_folder/spret_bt2_alignment
# processed_fastq_list=$trimmed_fastq_folder/sample_tn5_ligation.txt
# script=$script_folder/bowtie2_mapping.sh
# 
# echo "Start alignment using bowtie2..."
# echo input folder: $processed_fastq_folder
# echo processed_fastq_list: $processed_fastq_list
# echo output_folder: $bt2_output_folder
# mkdir -p $bt2_output_folder
# #start the alignment
# for sample in $(cat $processed_fastq_list); do echo Aligning $sample; echo "$script $sample $ref_bt2 $processed_fastq_folder/$sample.R1.trimmed.attached.fastq.gz $processed_fastq_folder/$sample.R2.trimmed.attached.fastq.gz $bt2_output_folder" | qsub -cwd -V -l mfree=20G -N bt2; done
# echo "All bt2 alignment done."
# 
############# split .bam file, use cut_off to define single cells
############# change bam_folder and split_bam_folder to bt2, mem_snp, bt2_snp and re-run
# echo
# echo "Start splitting .bam..."
# echo $(date)
# bam_folder=$all_output_folder/bwa_alignment
# split_bam_folder=$all_output_folder/split_bam_bwa
# script=$script_folder/split_bam.py
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# bam_sample_list=$trimmed_fastq_folder/sample_tn5_ligation.txt
# 
# echo input folder: $bam_folder
# echo output folder: $split_bam_folder
# echo round 1 barcode: $barcodes_tn5
# echo round 2 barcode: $barcodes_ligation
# echo cutoff: $cutoff
# 
# mkdir -p $split_bam_folder
# cat $bam_sample_list | while read line
# do
# mkdir -p $split_bam_folder/$line
# samtools view -H $bam_folder/$line.srt.bam > $split_bam_folder/$line/$line.header
# done
# 
# echo $cutoff
# #for sample in $(head -11 $bam_sample_list | tail -1); do echo Now splitting $sample; sem -j $core $python_path/python $script $bam_folder/$sample.srt.bam $split_bam_folder/$sample $sample $barcodes_tn5 $barcodes_ligation $cutoff; done
# for sample in $(cat $bam_sample_list); do echo Now splitting $sample; sem -j $core $python_path/python $script $bam_folder/$sample.srt.bam $split_bam_folder/$sample $sample $barcodes_tn5 $barcodes_ligation $cutoff; done
# sem --semaphoretimeout 1800
# echo
# echo "All sam file splitted."
# 
# mkdir -p $split_bam_folder/genome_coverage

#############################write summaries########################
bam_folder=$all_output_folder/bwa_alignment
split_bam_folder=$all_output_folder/split_bam_bwa
trimmed_fastq_folder=$all_output_folder/trimmed_fastq
bam_sample_list=$trimmed_fastq_folder/sample_tn5_ligation.txt

mkdir -p $split_bam_folder/aneuploidy
mkdir -p $split_bam_folder/aneuploidy/zero
mkdir -p $split_bam_folder/genome_coverage

for sample in $(cat $bam_sample_list)
do
# mkdir -p $split_bam_folder/$sample/bed_files
# rm $split_bam_folder/$sample/bed_files/* 
cat $split_bam_folder/$sample/$sample.sample_list.txt | while read BarComb
do
# $SAMTOOLS index $split_bam_folder/$sample/$BarComb.bam
# $SAMTOOLS idxstats $split_bam_folder/$sample/$BarComb.bam > $split_bam_folder/genome_coverage/$BarComb.idx
# echo $BarComb >> $split_bam_folder/genome_coverage/$sample.genome_coverage.txt
# awk '{if ($1 ~ /^chr/) print}' $split_bam_folder/genome_coverage/$BarComb.idx | awk '{sum += $3} END {print sum}' >> $split_bam_folder/genome_coverage/$sample.genome_coverage.txt
# 
awk 'FNR==NR{sum+=$3;next}; {print $1, $2, $3,sum}' $split_bam_folder/genome_coverage/$BarComb.idx{,} | grep -v "Un" | grep -v "random" | awk '{if($1 ~ /^chr/) print $1, $2, $3, $4}' > $split_bam_folder/aneuploidy/$BarComb.aneuploidy.txt
paste $split_bam_folder/aneuploidy/$BarComb.aneuploidy.txt $script_folder/mm10_chr_prop.txt | awk '{print $1, $2, $3, $3/$4, $4, ($3/$4)/$5}' > $split_bam_folder/aneuploidy/$BarComb.chr_prop.txt
awk '{if($6 < 0.1) print}' $split_bam_folder/aneuploidy/$BarComb.chr_prop.txt > $split_bam_folder/aneuploidy/zero/$BarComb.zero.txt
# 
# bedtools bamtobed -i $split_bam_folder/$sample/$BarComb.bam > $split_bam_folder/$sample/bed_files/$BarComb.bed
# grep "/1" $split_bam_folder/$sample/bed_files/$BarComb.bed > $split_bam_folder/$sample/bed_files/$BarComb.R1.bed
# awk '{if ($1 ~ /^chr/) print}' $split_bam_folder/$sample/bed_files/$BarComb.R1.bed > $split_bam_folder/$sample/bed_files/$BarComb.R1.human.bed
# awk '{if ($5 > 0) print}' $split_bam_folder/$sample/bed_files/$BarComb.R1.human.bed > $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.bed
# wc -l < $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.bed >> $split_bam_folder/genome_coverage/$sample.genome_coverage.txt
# awk '!seen[$1, $2]++' $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.bed > $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq1.bed
# awk '!seen[$1, $3]++' $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.bed > $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq2.bed
# bedtools intersect -a $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq1.bed -b $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq2.bed -f 1 -F 1 > $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq.bed
# rm $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq1.bed
# rm $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq2.bed
# wc -l < $split_bam_folder/$sample/bed_files/$BarComb.R1.human.Qgt0.uniq.bed >> $split_bam_folder/genome_coverage/$sample.genome_coverage.txt
done
done
# 
# mkdir -p $split_bam_folder/genome_coverage/single_cells
# 
# for sample in $(cat $bam_sample_list)
# do
# awk 'NR%4{printf "%s ",$0;next;}1' $split_bam_folder/genome_coverage/$sample.genome_coverage.txt | awk '{print $1, $2, $3, $4, $3/$4}' > $split_bam_folder/genome_coverage/single_cells/$sample.single_cells_all.txt
# done
# cat $split_bam_folder/genome_coverage/single_cells/*.single_cells_all.txt > $split_bam_folder/genome_coverage/single_cells/single_cell_summary_all.txt
# 
#############################generate vcf single cell: per chromosome#######################
# bulk_bam_folder=$all_output_folder/bwa_alignment/
# single_cell_bam_folder=$all_output_folder/split_bam_bwa/
# # bam_sample_list=$single_cell_bam_folder/all_single_cells/mouse.single_cells.bam_list.txt
# chr_file=$all_output_folder/chr.txt
# lib_file=$bulk_bam_folder/seq_lib.txt
# bulk_file=/net/shendure/vol8/projects/fly_recomb/nobackup/alignment/nextseq180206/sci_lianti_v2/output_lib190/bwa_alignment/mouse.2.srt.bam
# bulk_file_index=/net/shendure/vol8/projects/fly_recomb/nobackup/alignment/nextseq180206/sci_lianti_v2/output_lib190/bwa_alignment/mouse.2.srt.bam.bai
# output_folder=$bulk_bam_folder/vcf/
# 
# mkdir -p $output_folder
# mkdir -p $output_folder/bulk
# mkdir -p $output_folder/single_cells
# 
# cat $lib_file | while read -r lib num
# do
# mkdir -p $output_folder/single_cells/$lib
# cd $single_cell_bam_folder/$lib/
# ln -s $bulk_file mouse.srt.bam
# ln -s $bulk_file_index mouse.srt.bam.bai
# cat $chr_file | while read chr
# do
# echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/mm10.fa -C -q20 -Q20 -s2 -P20 -r $chr -L $num *.bam | gzip > $output_folder/single_cells/$lib/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=10G -N pileup
# # echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/dm3.fa -C -q20 -Q20 -s2 -P20 -r $chr -L $num $bulk_file $single_cell_bam_folder/$lib/*.bam | gzip > $output_folder/single_cells/$lib/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=2.5G -N pileup
# #echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m.fa -C -q20 -Q20 -s2 -P20 -r $chr -L0 $bulk_file | gzip > $output_folder/bulk/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=10G -N pileup
# #echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m.fa -C -q20 -Q20 -s1 -P20 -r $chr -L91 $bulk_file $single_cell_bam_folder/yi139*.bam | gzip > $output_folder/single_cells_cutoff1/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=5G -N pileup
# done
# done
# # 

#############################generate vcf single cell: per cell#######################
#############################get all variants in bulk, filter for absence in single cell######################
#############################bi-allelic filter#######################
# bulk_bam_folder=$all_output_folder/bwa_alignment/
# single_cell_bam_folder=$all_output_folder/split_bam_bwa/
# # bam_sample_list=$single_cell_bam_folder/all_single_cells/mouse.single_cells.bam_list.txt
# # chr_file=$all_output_folder/chr.txt
# lib_file=$bulk_bam_folder/seq_lib.txt
# bulk_file=/net/shendure/vol8/projects/fly_recomb/nobackup/alignment/nextseq180206/sci_lianti_v2/output_lib190/bwa_alignment/mouse.2.srt.bam
# bulk_file_index=/net/shendure/vol8/projects/fly_recomb/nobackup/alignment/nextseq180206/sci_lianti_v2/output_lib190/bwa_alignment/mouse.2.srt.bam.bai
# output_folder=$bulk_bam_folder/vcf/
# 
# mkdir -p $output_folder
# mkdir -p $output_folder/bulk
# mkdir -p $output_folder/single_cells
# mkdir -p $output_folder/by_cells
# 
# cat $lib_file | while read -r lib num
# do
# mkdir -p $output_folder/by_cells/$lib
# cd $single_cell_bam_folder/$lib/
# ln -s $bulk_file mouse.srt.bam
# ln -s $bulk_file_index mouse.srt.bam.bai
# cat $lib.sample_list.txt | while read sample
# do
# echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/mm10.fa -C -q1 -Q20 -s2 -P20 -L 1 mouse.srt.bam $sample.bam | grep -v "0,0:0,0:0$" | bcftools view -v snps -m 2 -M 2 -o $output_folder/by_cells/$lib/$sample.vcf.gz -O z" | qsub -cwd -V -l mfree=600M -N sec193
# # lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/mm10_natural_order.fa -C -q1 -Q20 -s2 -P20 -L 1 mouse.srt.bam $sample.bam | grep -v "0,0:0,0:0$" | bcftools view -v snps -m 2 -M 2 -o $output_folder/by_cells/$lib/$sample.vcf.gz -O z
# # echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/dm3.fa -C -q20 -Q20 -s2 -P20 -r $chr -L $num $bulk_file $single_cell_bam_folder/$lib/*.bam | gzip > $output_folder/single_cells/$lib/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=2.5G -N pileup
# #echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m.fa -C -q20 -Q20 -s2 -P20 -r $chr -L0 $bulk_file | gzip > $output_folder/bulk/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=10G -N pileup
# #echo "lianti-r109 pileup -cf /net/shendure/vol8/projects/fly_recomb/lh3/bwadb/hs37m.fa -C -q20 -Q20 -s1 -P20 -r $chr -L91 $bulk_file $single_cell_bam_folder/yi139*.bam | gzip > $output_folder/single_cells_cutoff1/raw.$chr.vcf.gz" | qsub -cwd -V -l mfree=5G -N pileup
# done
# done

#############################annotate variants with ref vcf##########################
# bulk_bam_folder=$all_output_folder/bwa_alignment/
# single_cell_bam_folder=$all_output_folder/split_bam_bwa/
# vcf_folder=$bulk_bam_folder/vcf/by_cells/
# annotated_vcf_folder=$bulk_bam_folder/vcf/by_cells/annotated
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# bam_sample_list=$trimmed_fastq_folder/sample_tn5_ligation.txt
# 
# cd $vcf_folder
# mkdir -p $annotated_vcf_folder
# 
# cat $bam_sample_list | while read lib
# do
# mkdir -p $annotated_vcf_folder/$lib
# cat $single_cell_bam_folder/$lib/$lib.sample_list.txt | while read cell
# do
# echo "java -jar ~/src/snpEff/SnpSift.jar annotate -id -noInfo  /net/shendure/vol8/projects/fly_recomb/lh3/vcf/SPRET_EiJ.v5.snps.dbSNP142.ucsc.includeLOW.vcf.gz $vcf_folder/$lib/$cell.vcf.gz | bgzip > $annotated_vcf_folder/$lib/$cell.filtered.vcf.gz" | qsub -cwd -V -l mfree=13G -N re193
# done
# done

############################test###############################################
# bulk_bam_folder=$all_output_folder/bwa_alignment/
# vcf_folder=$bulk_bam_folder/vcf/by_cells/
# script=$script_folder/sci-lianti/run-directory.R
# trimmed_fastq_folder=$all_output_folder/trimmed_fastq
# bam_sample_list=$trimmed_fastq_folder/sample_tn5_ligation_2.txt
# 
# cd $vcf_folder
# mkdir -p $vcf_folder/output2
# 
# cat $bam_sample_list | while read lib
# do
# echo $lib
# echo "Rscript --vanilla $script -i annotated/$lib/ -or output2/$lib.rds -hb output2/$lib.hb.tab -mb output2/$lib.mb.tab -hp output2/$lib.hp.pdf -mp output2/$lib.mp.pdf -c output2/$lib.cellstatus.tab" | qsub -cwd -V -l mfree=27G -N R193
# # echo "Rscript --vanilla $script $lib/ output/$lib.rds output/$lib.bed output/$lib.pdf output/$lib.tab" | qsub -cwd -V -l mfree=3G -N bt186
# done
