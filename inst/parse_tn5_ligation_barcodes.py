import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial


'''
    this script accepts and opens trimmed R1 and ordered R2 files splited according to sss barcodes, an reformatted info file from cutadapt, a bridge sequence, a tn5 barcode list and a ligation barcode list; 
    then extracts and matches the nearest tn5 and ligation barcodes and attach barcodes to R1 and R2 read names.
    the output files include new R1 and R2 files, noME.R1.fq.gz
    
'''

def tn5_ligation_barcode_attach_list(sample, input_folder, barcode_tn5, barcode_ligation, mismatch_tn5, mismatch_ligation):
    Read1 = input_folder + "/" + sample + ".R1.fq.gz"
    Read2 = input_folder + "/" + sample + ".R2.fq.gz"
    R1_out_file = input_folder + "/" + sample + ".R1.trimmed.attached.fastq.gz"
    R2_out_file = input_folder + "/" + sample + ".R2.attached.fastq.gz"
    bc1_bc2_file = input_folder + "/" + sample + ".R1.bc1.bc2.txt.gz"
    noME_R1_out_file = input_folder + "/bkp/" + sample + ".noME.R1.fq.gz"
    mismatch_tn5 = int(mismatch_tn5)
    mismatch_ligation = int(mismatch_ligation)
    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    f3 = gzip.open(bc1_bc2_file)
    f4 = gzip.open(R1_out_file, 'wb')
    f5 = gzip.open(R2_out_file, 'wb')
    f6 = gzip.open(noME_R1_out_file, 'wb')
    
    line1 = f1.readline()
    line2 = f2.readline()
    line3 = f3.readline()
    
    while (line1):
        name_r1=line1
        name_r2=line2
        barcode_target = line3.split()
        if len(barcode_target) < 2:
            barcode_target = ["NA","NA"]
        target_tn5 = barcode_target[1]
        target_ligation = barcode_target[0]
        find = False
        if (target_tn5 != "NA"):
            for tn5 in barcode_tn5:
                mm_tn5 = distance(tn5, target_tn5)
                if (mm_tn5 <= mismatch_tn5):
                    for ligation in barcode_ligation:
                        mm_ligation = distance(ligation, target_ligation)
                        if(mm_ligation <= mismatch_ligation):
                            find = True
                            first_line_r1 = '@' + tn5 + ',' + ligation + ',' + name_r1[1:]
                            first_line_r2 = '@' + tn5 + ',' + ligation + ',' + name_r2[1:]
                            f4.write(first_line_r1)
                            f5.write(first_line_r2)
                            
                            second_line_r1 = f1.readline()
                            second_line_r2 = f2.readline()
                            f4.write(second_line_r1)
                            f5.write(second_line_r2)
                            
                            third_line_r1 = f1.readline()
                            third_line_r2 = f2.readline()
                            f4.write(third_line_r1)
                            f5.write(third_line_r2)
                            
                            four_line_r1 = f1.readline()
                            four_line_r2 = f2.readline()
                            f4.write(four_line_r1)
                            f5.write(four_line_r2)
                            
                            line1 = f1.readline()
                            line2 = f2.readline()
                            line3 = f3.readline()
                            break
                    if (find == True):
                        break
            if (find == False):
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line3 = f3.readline()
        #write noME file
        else: 
            first_line_r1 = name_r1
            first_line_r2 = name_r2
            f6.write(first_line_r1)
            second_line_r1 = f1.readline()
            second_line_r2 = f2.readline()
            f6.write(second_line_r1)
            third_line_r1 = f1.readline()
            third_line_r2 = f2.readline()
            f6.write(third_line_r1)
            four_line_r1 = f1.readline()
            four_line_r2 = f2.readline()
            f6.write(four_line_r1)
            line1 = f1.readline()
            line2 = f2.readline()
            line3 = f3.readline()
            
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()

# this function accept an input folder and a output folder and then generate the output file with the index
def attach_tn5_ligation_barcode_files(input_folder, sample_file, barcode_tn5_file, barcode_ligation_file, core_number, mismatch_tn5, mismatch_ligation):
    
    init_message = '''
    --------------------------start attaching barcodes-----------------------------
    input and output folder: %s
    core number: %s
    sample file: %s
    barcode files: %s, %s
    ___________________________________________________________________________
    ''' %(input_folder, core_number, sample_file, barcode_tn5_file, barcode_ligation_file)
    
    print(init_message)
    
    # generate the barcode list:
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
    
    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_tn5 and barcode_ligation to run tn5_ligation_barcode_attach_list
    sampleID = open(sample_file)
    sample_list = []
    for line in sampleID:
        sample = line.strip()
        sample_list.append(sample)
    sampleID.close()
   
#    tn5_ligation_barcode_attach_list(sample_list[0], input_folder, barcode_tn5, barcode_ligation, mismatch_tn5, mismatch_ligation) 

    # parallele for the functions
    p = Pool(processes = int(core_number))
    #print("Processing core number: ", core_number)
    func = partial(tn5_ligation_barcode_attach_list, input_folder=input_folder, barcode_tn5=barcode_tn5, barcode_ligation=barcode_ligation, mismatch_tn5=mismatch_tn5, mismatch_ligation=mismatch_ligation)
    result = p.map(func, sample_list)
    p.close()
    p.join()

    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~barcode attachment done~~~~~~~~~~~~~~~~~~'''
    print com_message
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sample_file = sys.argv[2]
    barcode_tn5_file = sys.argv[3]
    barcode_ligation_file = sys.argv[4]
    core = sys.argv[5]
    mismatch_tn5 = sys.argv[6]
    mismatch_ligation = sys.argv[7]
    attach_tn5_ligation_barcode_files(input_folder, sample_file, barcode_tn5_file, barcode_ligation_file, core, mismatch_tn5, mismatch_ligation)
