#####################################################################
#####################################################################

## (1) Combine four fastq.gz files with the same header, but four different flow cell lane numbers
## into single fastq.gz file. This works for all fastq.gz files in a designated folder

## Go to scripts
## python <script> <folder that contains all fastq.gz files>
## python GetSingleReadFastq.py ../fastq/Batch_retest/

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
parser.add_argument('output_folder', help='Folder name for output files',type=str)

# Saves input to args object
args = parser.parse_args()
                
#####################################################################
##################### Get names in the folder #######################
#####################################################################
file_out = args.output_folder
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
    if 'Undetermined' != header:        
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
        
print('Finished')