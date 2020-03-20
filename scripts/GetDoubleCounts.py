#####################################################################
#####################################################################

## Go to /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/scripts/
## python <script> <folder that contains all fastq.gz files> <output folder>
## python GetDoubleCounts.py /mnt/lab_data/bassik/kyuhohan/NextSeq/bcl2fastq/KRAS_PPI_20by20_180206_Y_Y_I/ /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/counts/20180214_KRAS_PPI_20by20/ KRAS_20by20 KRAS_20by20
## python GetDoubleCounts.py /mnt/lab_data/bassik/kyuhohan/NextSeq/bcl2fastq/AmyStewart_Y_Y_I/ /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/counts/AmyStewart_20180227/ amy_short_mu6 amy_long_hu6
## python GetDoubleCounts.py /mnt/lab_data/bassik/kyuhohan/NextSeq/bcl2fastq/KRAS_pgRNA_Kaitlyn_180316/ /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/counts/20180316_KRAS_PPI_pgRNA/ KRAS_PPI_pgRNA KRAS_PPI_pgRNA

#####################################################################
#####################################################################

import subprocess
import shlex
import sys
import csv
import os
import argparse
import HTSeq
from collections import defaultdict
import time

def sorteddict2tab(dict,file):
    tableWriter = csv.writer(open(file, 'w'), delimiter='\t')
    for key in sorted(dict):
        tableWriter.writerow([key, dict[key]])

def sorteddict2tabline(dict,file):
    tableWriter = csv.writer(open(file, 'w'), delimiter='\t')
    for key in sorted(dict):
        tableWriter.writerow([key] + dict[key])

def execute(command):
    print command
    subprocess.call(shlex.split(command))

# Initiates argument parser
parser = argparse.ArgumentParser(description='Align FASTQ and make count file')

# Non-optional arguments:
parser.add_argument('fastq_folder', help='Folder path for input fastq files',type=str)
parser.add_argument('count_folder', help='Folder name for output count files',type=str)
parser.add_argument('index1_name', help='index1_name', type=str)
parser.add_argument('index2_name', help='index2_name', type=str)

# optional arguments
parser.add_argument('-m', '--mismatch',help='No. of mismatches allowed (1-3). Default is 1', type=int, default=1)

# Saves input to args object
args = parser.parse_args()
                
#####################################################################
##################### Get names in the folder #######################
#####################################################################
file_out = args.count_folder
seq_folder, seq_name = os.path.split(args.fastq_folder)

file_headers = []

for root, dirs, files in os.walk(seq_folder):
    for fil in files:
        if fil.endswith('fastq.gz'):            
            Temp = fil.split('_')[:-4]
            header = ''
            for tmp in Temp:
                header=header+tmp
                header=header+'_'
            header = header[:-1]
            
            file_headers.append(header)
        
file_headers = list(set(file_headers))

#####################################################################
### concaternate, unzip, align gz files,  ###########################
### and make counts with the same headers ###########################
#####################################################################
bowtie='/software/bowtie/1.1.1/bowtie'

for header in file_headers:
    R1_files=[]
    R2_files=[]
    for root, dirs, files in os.walk(seq_folder):
        for fil in files: 
            if fil.startswith(header):
                if 'R1' in fil[5:]:
                    R1_files.append(os.path.join(root, fil))
        
                elif 'R2' in fil[5:]:
                    R2_files.append(os.path.join(root, fil))
        
                elif 'I1' in fil[5:] or 'I2' in fil[5:]:
                    pass
        
                else:
                    print('Warning: conflicting file names:\n' + fil)
        
    print('Read 1 files:\n' + '\n'.join(R1_files))
    print('Read 2 files:\n' + '\n'.join(R2_files))
    
    # 1) concaternate
    try:
        subprocess.check_call('cat ' + ' '.join(R1_files) + ' > '
                            + file_out + header + '_R1_all.fastq.gz', shell=True)
    except:
        sys.exit('Shell error')

    try:
        subprocess.check_call('cat ' + ' '.join(R2_files) + ' > '
                                + file_out + header + '_R2_all.fastq.gz', shell=True)
    except:
        sys.exit('Shell error')
        
    print('Unzipping read one'+' : '+header)
    
    # 2) unzip
    try:
        subprocess.check_call('gunzip ' + file_out + header + '_R1_all.fastq.gz -f', shell=True)
    except:
        sys.exit('Shell error')
        
    try:
        subprocess.check_call('gunzip ' + file_out + header + '_R2_all.fastq.gz -f', shell=True)
    except:
        sys.exit('Shell error')

    # 3) bowtie alignment        
    fastq1 = file_out + header + '_R1_all.fastq'
    fastq2 = file_out + header + '_R2_all.fastq'    
    
    execute(bowtie+' -v ' + str(args.mismatch) + ' -a -p 8 --trim3 3 /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/indices/'+args.index1_name+' -q '+ fastq1 + ' ' + file_out + header + '_r1.map --un ' + file_out + header + '_r1.unmapped')  
    execute(bowtie+' -v ' + str(args.mismatch) + ' -a -p 8 --trim3 3 /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/indices/'+args.index2_name+' -q '+ fastq2 + ' ' + file_out + header + '_r2.map --un ' + file_out + header + '_r2.unmapped')  
    
    # 4) make counts
    file1=open(file_out + header + '_r1.map','rb')
    file2=open(file_out + header + '_r2.map','rb')

    fileread_r1 = csv.reader(file1, delimiter='\t')
    fileread_r2 = csv.reader(file2, delimiter='\t')
    
    # save unique identifier of each paired end read and its corresponding shRNA in dictionary 
    R1UniqId2shRNA={}
    R2UniqId2shRNA={}
    TotalReads1 = 0
    TotalReads2 = 0
    R1count={}
    R2count={}
    
    for line in fileread_r1:
        TotalReads1+=1
        if line[1]=='+':
            
            temp_uniq=line[0].split(' ')[0]
            if temp_uniq in R1UniqId2shRNA:
                continue
            else:
                R1UniqId2shRNA[temp_uniq]=line[2]
            if line[2] in R1count:
                R1count[line[2]]+=1
            else:
                R1count[line[2]]=1
                        
        if TotalReads1%1000000 == 0:
            print 'Reads analyzed so far: '+str(TotalReads1/1000000)+' million'
        
    print 'Total number of reads mapped to '+sys.argv[3]+': ', str(TotalReads1)
    file1.close()
    
    for line in fileread_r2:
        TotalReads2+=1
        if line[1]=='+':
            temp_uniq=line[0].split(' ')[0]
            if temp_uniq in R2UniqId2shRNA:
                continue
            else:
                R2UniqId2shRNA[temp_uniq]=line[2]
            if line[2] in R2count:
                R2count[line[2]]+=1
            else:
                R2count[line[2]]=1
    
        if TotalReads2%1000000 == 0:
            print 'Reads analyzed so far: '+str(TotalReads2/1000000)+' million'
            
    print 'Total number of reads mapped to '+sys.argv[3]+': ', str(TotalReads2)
    file2.close()
    
    double_shRNA_count={}
    TotalReads=0
    for key, value in R1UniqId2shRNA.iteritems():
        TotalReads+=1
        if key in R2UniqId2shRNA:
            temp_double=value+'____'+R2UniqId2shRNA[key]
            if temp_double in double_shRNA_count:
                double_shRNA_count[temp_double]+=1
            else:
                double_shRNA_count[temp_double]=1                        
        if TotalReads%1000000 == 0:
            print 'Reads analyzed so far: '+str(TotalReads/1000000)+' million'                        
    print 'Total number of unique double elements : ', str(len(double_shRNA_count))
                        
    sorteddict2tab(R1count,file_out+header+'_r1.counts')
    sorteddict2tab(R2count,file_out+header+'_r2.counts')
    sorteddict2tab(double_shRNA_count,file_out+header+'_double.counts')
    
    # 5) Compresses and/or erases files using linux shell commands

    print('Deleting gathered reads')
    try:
        subprocess.check_call('rm ' + file_out + header + '_R1_all.fastq', shell=True)
        subprocess.check_call('rm ' + file_out + header + '_R2_all.fastq', shell=True)
    except:
        sys.exit('Shell error')
    
    print('Compressing mapped reads')
    try:
        subprocess.check_call('gzip ' + file_out + header + '_r1.map -f', shell=True)
        subprocess.check_call('gzip ' + file_out + header + '_r2.map -f', shell=True)
    except:
        sys.exit('Shell error')
    
    print('Compressing unmapped reads')
    try:
        subprocess.check_call('gzip ' + file_out + header + '_r1.unmapped -f', shell=True)
        subprocess.check_call('gzip ' + file_out + header + '_r2.unmapped -f', shell=True)
    except:
        sys.exit('Shell error')
    
print('Finished')