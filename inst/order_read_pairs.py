import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial


'''
    this script accepts and opens a R1 file, a R2 file, an RT primer sequence, and a SSS barcode list; then reorders the read pairs, and extracts and matches the nearest sss barcode and attach the sss barcode to R1 and R2 read names.
    the output files include ordered R1 and R2 files, noRTprimer.R1.fq.gz and noRTprimer.R2.fq.gz
    
'''

def SSS_barcode_attach_list(sample, input_folder, output_folder, barcode_list, RT_primer, mismatch_sss, mismatch_RT):
    #open the read1, read2, and output files
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R2.fastq.gz"
    R1_out_file = output_folder + "/" + sample + ".R1.ordered.fastq.gz"
    R2_out_file = output_folder + "/" + sample + ".R2.ordered.fastq.gz"
    noRT_R1_out_file = output_folder + "/" + sample + ".noRTprimer.R1.fq.gz"
    noRT_R2_out_file = output_folder + "/" + sample + ".noRTprimer.R2.fq.gz"
    mismatch_sss = int(mismatch_sss)
    mismatch_RT = int(mismatch_RT)
    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    f3 = gzip.open(R1_out_file, 'wb')
    f4 = gzip.open(R2_out_file, 'wb')
    f5 = gzip.open(noRT_R1_out_file, 'wb')
    f6 = gzip.open(noRT_R2_out_file, 'wb')
    
    line1 = f1.readline()
    line2 = f2.readline()

    while (line1):
        name_r1=line1
        name_r2=line2
        line1 = f1.readline()
        line2 = f2.readline()
        target_rt_r1 = line1[9:29]
        target_rt_r2 = line2[9:29]
        target_sss_r1 = line1[4:10]
        target_sss_r2 = line2[4:10]
        
        junk=True
        dist_rt_r1 = distance(RT_primer, target_rt_r1)
        dist_rt_r2 = distance(RT_primer, target_rt_r2)
        
        if (dist_rt_r1 <= mismatch_RT):
            junk=False
            for barcode in barcode_list:
                mismatch = distance(barcode, target_sss_r1)
                find = False
                
                if (mismatch <= mismatch_sss):
                    find = True
                    UMI = line1[:4]
                    first_line_r1 = '@' + barcode + ',' + UMI + ',' + name_r1[1:]
                    first_line_r2 = '@' + barcode + ',' + UMI + ',' + name_r2[1:]
                    f3.write(first_line_r1)
                    f4.write(first_line_r2)
    
                    second_line_r1 = line1
                    second_line_r2 = line2
                    f3.write(second_line_r1)
                    f4.write(second_line_r2)
    
                    third_line_r1 = f1.readline()
                    third_line_r2 = f2.readline()
                    f3.write(third_line_r1)
                    f4.write(third_line_r2)
    
                    four_line_r1 = f1.readline()
                    four_line_r2 = f2.readline()
                    f3.write(four_line_r1)
                    f4.write(four_line_r2)
    
                    line1 = f1.readline()
                    line2 = f2.readline()
                    break
                    
            if find == False:
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                
        #swap R1 R2 but keep read name
        elif (dist_rt_r2 <= mismatch_RT):
            junk=False
            for barcode in barcode_list:
                mismatch = distance(barcode, target_sss_r2)
                find = False
                
                if (mismatch <= mismatch_sss):
                    find = True
                    UMI = line2[:4]
                    first_line_r1 = '@' + barcode + ',' + UMI + ',' + 'swap' + ',' + name_r1[1:]
                    first_line_r2 = '@' + barcode + ',' + UMI + ',' + 'swap' + ',' + name_r2[1:]
                    f3.write(first_line_r1)
                    f4.write(first_line_r2)
    
                    second_line_r1 = line2
                    second_line_r2 = line1
                    f3.write(second_line_r1)
                    f4.write(second_line_r2)
    
                    third_line_r1 = f2.readline()
                    third_line_r2 = f1.readline()
                    f3.write(third_line_r1)
                    f4.write(third_line_r2)
    
                    four_line_r1 = f2.readline()
                    four_line_r2 = f1.readline()
                    f3.write(four_line_r1)
                    f4.write(four_line_r2)
    
                    line1 = f1.readline()
                    line2 = f2.readline()
                    break
                    
            if find == False:
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                
        #write junk file
        else: 
            junk = True
            first_line_r1 = name_r1
            first_line_r2 = name_r2
            f5.write(first_line_r1)
            f6.write(first_line_r2)

            second_line_r1 = line1
            second_line_r2 = line2
            f5.write(second_line_r1)
            f6.write(second_line_r2)

            third_line_r1 = f1.readline()
            third_line_r2 = f2.readline()
            f5.write(third_line_r1)
            f6.write(third_line_r2)

            four_line_r1 = f1.readline()
            four_line_r2 = f2.readline()
            f5.write(four_line_r1)
            f6.write(four_line_r2)
            line1 = f1.readline()
            line2 = f2.readline()
            
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()

# this function accept an input folder and a output folder and then generate the output file with the index
def attach_SSS_barcode_files(input_folder, sampleID, output_folder, barcode_sss, core_number, mismatch_sss, mismatch_RT, RT_primer):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    barcode file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, barcode_sss)
    
    print(init_message)
    
    # generate the barcode list:
    barcode_list = []
    barcodes = open(barcode_sss)
    for barcode in barcodes:
        barcode_list.append(barcode.strip())
    barcodes.close()
    
    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run SSS_barcode_attach_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core_number))
    #print("Processing core number: ", core_number)
    func = partial(SSS_barcode_attach_list, input_folder = input_folder, output_folder=output_folder, barcode_list=barcode_list, mismatch_sss = mismatch_sss, mismatch_RT=mismatch_RT, RT_primer=RT_primer)

    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~SSS and UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print com_message
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    barcode_sss = sys.argv[4]
    core = sys.argv[5]
    mismatch_sss = sys.argv[6]
    mismatch_RT = sys.argv[7]
    RT_primer = sys.argv[8]
    attach_SSS_barcode_files(input_folder, sampleID, output_folder, barcode_sss, core, mismatch_sss, mismatch_RT, RT_primer)