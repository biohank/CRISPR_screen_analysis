#####################################################################
#####################################################################

## (2) Align fastq.gz files against an index file of a sgRNA library to generate count files

## Go to scripts
## python GetSingleCounts.py <folder with fastq.gz files processed in step (1)> <output file> <header of bowtie index file in "indices" folder>
## ex) python GetSingleCounts.py ../fastq/Batch_retest/ ../counts/Batch_retest/ Lung3D_Retest

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
                elif 'I1' in fil[5:] or 'I2' in fil[5:]:
                    pass
        
                else:
                    print('Warning: conflicting file names:\n' + fil)
        
    print('Read 1 files:\n' + '\n'.join(R1_files))
    
    # 1) concaternate
    try:
        subprocess.check_call('cat ' + ' '.join(R1_files) + ' > '
                            + file_out + header + '_R1_all.fastq.gz', shell=True)
    except:
        sys.exit('Shell error')
        
    print('Unzipping read one'+' : '+header)
    
    # 2) unzip
    try:
        subprocess.check_call('gunzip ' + file_out + header + '_R1_all.fastq.gz -f', shell=True)
    except:
        sys.exit('Shell error')

    # 3) bowtie alignment        
    fastq1 = file_out + header + '_R1_all.fastq'
    
    execute(bowtie+' -v 1 -a -p 8 --trim3 3 /mnt/lab_data/bassik/kyuhohan/Kyuho_Screening_Analysis/indices/'+args.index1_name+' -q '+ fastq1 + ' ' + file_out + header + '_r1.map --un ' + file_out + header + '_r1.unmapped')  
    
    # 4) make counts
    file1=open(file_out + header + '_r1.map','rb')

    fileread_r1 = csv.reader(file1, delimiter='\t')
    
    # save unique identifier of each paired end read and its corresponding shRNA in dictionary 
    R1UniqId2shRNA={}
    TotalReads1 = 0
    R1count={}
    
    for line in fileread_r1:
        TotalReads1+=1
        if line[1]=='-':
            
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
                        
    sorteddict2tab(R1count,file_out+header+'_r1.counts')
    
    # 5) Compresses and/or erases files using linux shell commands

    print('Deleting gathered reads')
    try:
        subprocess.check_call('rm ' + file_out + header + '_R1_all.fastq', shell=True)
    except:
        sys.exit('Shell error')
    
    print('Compressing mapped reads')
    try:
        subprocess.check_call('gzip ' + file_out + header + '_r1.map -f', shell=True)
    except:
        sys.exit('Shell error')
    
    print('Compressing unmapped reads')
    try:
        subprocess.check_call('gzip ' + file_out + header + '_r1.unmapped -f', shell=True)
    except:
        sys.exit('Shell error')
    
print('Finished')