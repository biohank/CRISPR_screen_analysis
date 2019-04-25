#####################################################################################
# Kyuho Han (Some codes are adapted from David Morgen's scripts)
# 10/10/2017
# Compare two count files to calculate normalized log fold enrichment of sgRNAs
#####################################################################################

#run GetGIMatrix.py ../counts/L_Plas_r1.counts ../counts/L_Plas_r2.counts ../results/1521_4L_t0_vs_invivo_combo ../results/1521_4L_t0_vs_invivo_combo

from __future__ import division

# import matplotlib as mpl 
# mpl.use('Agg')
# mpl.rcParams['pdf.fonttype'] = 42
# import matplotlib.pyplot as plt
# plt.close("all")
# import matplotlib.cm as cm

import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as st
import pickle        

from scipy.interpolate import interp1d

import subprocess
import shlex
import sys
import csv
import os
import argparse
from FunctionLib import *

#####################################################################################
# Version number
current_version = '0.1'

#####################################################################################
# Parse input
parser = argparse.ArgumentParser(description='Make matrix files for GI map')

# Non-optional arguments
parser.add_argument('read1', help='File for read1 count', type=str)
parser.add_argument('read2', help='File for read2 counts', type=str)
parser.add_argument('input', help='header for input file', type=str)
parser.add_argument('output', help='Name for output file', type=str)

# Options for element IDs
parser.add_argument('-fs', '--force_symmetric',help='Forcing GI map to be symmetric. Default is n', type=str, default='n', choices=['y','n'])
parser.add_argument('-nc', '--neg_ctrl',help='Symbol used to denote negative controls. Default is safe.',type=str, default='safe')
parser.add_argument('-rm', '--remove_ctrl',help='Remove negative controls from GI matrix. Default is y',type=str, default='y', choices=['y','n'])
parser.add_argument('-sn', '--split_single',help='Delimiter for parsing single element name. Default is _',type=str, default='_')
parser.add_argument('-gi', '--gene_idx',help='Position of genename in parsed element name. Default is 1', type=int, default=1)

# Options for computation
parser.add_argument('-p', '--proccessors', dest='nums',help='Number of proccessors to use; default is 20', type=int,default=20)

args = parser.parse_args()

###############################################################################
# Processes and checks input arguments

## Check if you have read1 file

file_in = args.read1
try:
    with open(file_in, 'r') as in_open:
        print('read1 conut file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter='\t')
        front_genes = []
        for line in filein_csv:
            tmp_gene = line[0].split(args.split_single)[args.gene_idx]
            front_genes.append(tmp_gene)
        front_genes = list(set(front_genes))
except:
    sys.exit('Cannot find the following count file:\n' + file_in)

## Check if you have read2 file

file_in = args.read2
try:
    with open(file_in, 'r') as in_open:
        print('read2 conut file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter='\t')
        rear_genes = []
        for line in filein_csv:
            tmp_gene = line[0].split(args.split_single)[args.gene_idx]
            rear_genes.append(tmp_gene)
        rear_genes = list(set(rear_genes))
except:
    sys.exit('Cannot find the following count file:\n' + file_in)
    
## Check if you have input phenotype files

file_in = args.input+'_phenotype.csv'
try:
    with open(file_in, 'r') as in_open:
        print('input phenotype file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter=',')
        
        gitscore_dic = {}
        gimscore_dic = {}
        giqscore_dic = {}
        
        for i,line in enumerate(filein_csv):
            if i==0:
                continue
                
            tmp_gp = line[0]
            tmp_tscore = float(line[13])
            tmp_mscore = -1 * np.log10(float(line[15])) * np.sign(tmp_tscore)
            tmp_qscore = -1 * np.log10(float(line[16])) * np.sign(tmp_tscore)

            gitscore_dic[tmp_gp]=tmp_tscore
            gimscore_dic[tmp_gp]=tmp_mscore
            giqscore_dic[tmp_gp]=tmp_qscore           
                       
except:
    sys.exit('Cannot find the following phenotype file:\n' + file_in)

## Generate matrix and export matrix files

if args.force_symmetric == 'y':
    all_genes = front_genes + rear_genes
    all_genes = list(set(all_genes))
    
    front_genes = []
    rear_genes = []
    for g in all_genes:
        front_genes.append(g)
        rear_genes.append(g)
    
if args.remove_ctrl == 'y':
    front_genes.remove(args.neg_ctrl)
    rear_genes.remove(args.neg_ctrl)

tscore_mat = pd.DataFrame(data=None,index=front_genes,columns=rear_genes)
mscore_mat = pd.DataFrame(data=None,index=front_genes,columns=rear_genes)
qscore_mat = pd.DataFrame(data=None,index=front_genes,columns=rear_genes)
    
for i,fg in enumerate(front_genes):
    for j,rg in enumerate(rear_genes):
        name_front = np.sort([fg,rg])[0]
        name_rear = np.sort([fg,rg])[1]
        name_gp = name_front+'__'+name_rear
        
        if name_gp in gitscore_dic:
            tval = gitscore_dic[name_gp]
        else:
            tval = NaN

        if name_gp in gimscore_dic:
            mval = gimscore_dic[name_gp]
        else:
            mval = NaN
                        
        if name_gp in giqscore_dic:
            qval = giqscore_dic[name_gp]
        else:
            qval = NaN
                
        tscore_mat.iloc[i,j]=tval
        mscore_mat.iloc[i,j]=mval
        qscore_mat.iloc[i,j]=qval

file_out = args.output

tscore_mat.to_csv(file_out+'_tscore_mat.txt',sep='\t')
mscore_mat.to_csv(file_out+'_mscore_mat.txt',sep='\t')
qscore_mat.to_csv(file_out+'_qscore_mat.txt',sep='\t')