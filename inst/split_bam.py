import subprocess
import sys
import pysam
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


'''
    this script accepts and opens trimmed R1 and ordered R2 files splited according to sss barcodes, an reformatted info file from cutadapt, a bridge sequence, a tn5 barcode list and a ligation barcode list; 
    then extracts and matches the nearest tn5 and ligation barcodes and attach barcodes to R1 and R2 read names.
    the output files include new R1 and R2 files, noME.R1.fq.gz
    
'''

def bam_barcode_count(bam_file, barcode_tn5_file, barcode_ligation_file):
    
    #generate the barcode list and barcode dictionary
    barcode_tn5 = []
    barcodes = open(barcode_tn5_file)
    for barcode in barcodes:
        barcode_tn5.append(barcode.strip())
    barcodes.close()
    
    barcode_ligation = []
    barcodes = open(barcode_ligation_file)
    for barcode in barcodes:
        barcode_ligation.append(barcode.strip())
    barcodes.close()
    
    barcode_combination = []
    barcode_dic = {}
    for tn5 in barcode_tn5:
        for ligation in barcode_ligation:
            barcode = tn5 + ligation
            barcode_combination.append(barcode)
            barcode_dic[barcode] = 0
    #read the sam file, and count the number per barcode
    input_sam = pysam.AlignmentFile(bam_file, "rb")
    for line in input_sam.fetch():
#        if (line[0] == '@'):
        if (False):
            continue
        else:
            barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            barcode_dic[barcode] += 1
    input_sam.close()
    return barcode_dic

def split_bam(bam_file, output_folder, sample_file, barcode_tn5_file, barcode_ligation_file, cut_off):
    '''
    this script accept a bam file, a sample name, an output folder, two barcode files and a cut_off value,
    then it will call the bam_barcode_count function and get the total read count per barcode,
    then it use the cut_off value to filter the barcode,
    and generate the output samfile for single cells,
    generate the reads distribution in the split_bam/read_distribution_barcode;
    '''
    
    # generate the count per barcode
    barcode_count = bam_barcode_count(bam_file, barcode_tn5_file, barcode_ligation_file)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (bam_file.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('frequency')
    plt.xlabel('Number of unique reads')
    fig_output = output_folder + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    read_dist = open(output_folder + '/' + plot_name + '.txt', 'w')
    for barcode in barcode_count:
        line = barcode + ', %d\n' %(barcode_count[barcode])
        read_dist.write(line)
    read_dist.close()

    #Generate the read distribution in the split_bam/read_distribution_barcode
    
    #filter the barcode based on the cut_off value
    barcode_filtered = []
    for barcode in barcode_count:
        if (barcode_count[barcode] >= int(cut_off)):
            barcode_filtered.append(barcode)

    #generate the output sam file and sample_list file
    sample_list_file = open(output_folder + '/' + plot_name + '.' + 'sample_list.txt', 'w')
    output_files = {}
    paired_reads = {}
    for barcode in barcode_filtered:
        output_file = output_folder + '/' + plot_name + '.' + barcode + '.bam'
        output_files[barcode] = open(output_file, 'w')
        sample_list_file.write(plot_name + '.' + barcode + '\n')
    
    # output each read to the output sam file
    input_sam = pysam.AlignmentFile(bam_file, "rb")
    for barcode in barcode_filtered:
        outname=output_files[barcode]
        paired_reads[outname] = pysam.AlignmentFile(outname, "wb", template=input_sam)
        
    for line in input_sam.fetch():
#        if (line[0] == '@'):
        if(False):
            for barcode in barcode_filtered:
                output_files[barcode].write(line)
        else:
            barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            if barcode in barcode_filtered:
                outname=output_files[barcode]
                paired_reads[outname].write(line)
#                output_files[barcode].write(str(line))
    for barcode in barcode_filtered:
        outname=output_files[barcode]
        paired_reads[outname].close()
    
    #close the files:
    sample_list_file.close()
    input_sam.close()
    for barcode in barcode_filtered:
        output_files[barcode].close()

if __name__ == "__main__":
    bam_file = sys.argv[1]
    output_folder = sys.argv[2]
    sample_file = sys.argv[3]
    barcode_tn5_file = sys.argv[4]
    barcode_ligation_file = sys.argv[5]
    cut_off = sys.argv[6]
    split_bam(bam_file, output_folder, sample_file, barcode_tn5_file, barcode_ligation_file, cut_off)
