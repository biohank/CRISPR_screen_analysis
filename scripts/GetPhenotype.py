#####################################################################################
# Kyuho Han (Some codes are adapted from David Morgen's scripts)
# 10/10/2017
# Compare two count files to calculate normalized log fold enrichment of sgRNAs
#####################################################################################

# run GetPhenotype.py ../results/20171031_KRAS_3D/T0_vs_D19_2D_U1
# run GetPhenotype.py ../results/20171208_Youngtae/D16_CONT_A_vs_D16_IRA_A -m

# run GetPhenotype.py ../results/20171009_SL_Screens/Plas_vs_D15_Unt1
# run GetPhenotype.py ../results/20171009_SL_Screens/Plas_vs_D15_Unt2

from __future__ import division

import matplotlib as mpl 
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
plt.close("all")

import matplotlib.cm as cm

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
parser = argparse.ArgumentParser(description='Get phenotypes of genes (pZ)')

# Non-optional arguments
parser.add_argument('file_head', help='file head for enrichment and record input file and for output file', type=str)

# Options for analysis
parser.add_argument('-s', '--statistics',help='Additional statistics to use. Tscore will be calculated regardless of this option. Default is none', type=str, default='none', choices=['bootstrap','mwu','none','all'])
parser.add_argument('-br', '--bootstrap_repeat',help='No. of repeats for bootstrap. Default is 1000', type=int, default=1000)
parser.add_argument('-c', '--percentile_cut',help='Minimum percentile to calculate single ko. Default is 0.5', type=float, default=0.5) # not implemented yet
parser.add_argument('-f', '--fitting',help='data fitting to either all data or to neg_ctrl. Default is neg_ctrl', type=str, default='neg_ctrl', choices=['all','neg_ctrl'])
parser.add_argument('-g', '--gi_option',help='Sign change for GI depending on the exp growth phenotype. Default is corr_sign', type=str, default='corr_sign', choices=['corr_sign','unch_sign'])

# Options for additional info files
parser.add_argument('-m', '--mouse', action='store_true', help='Uses mouse gene information.')

# Options for computation
parser.add_argument('-p', '--proccessors', dest='nums',help='Number of proccessors to use; default is 20', type=int,default=20)
args = parser.parse_args()

###############################################################################
# Processes and checks input arguments

## Check if you have record files for enrichment files

file_in = args.file_head + '_record.txt'
try:
    with open(file_in, 'r') as in_open:
        print('Record file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter='\t')
        record_dic = {}
        for line in filein_csv:
            record_dic[line[0]]=line[1]   
except:
    sys.exit('Cannot find the following record file:\n' + file_in)

file_in = args.file_head + '_enrichment.csv'
try:
    with open(file_in, 'r') as in_open:
        print('Enrichment file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter=',')
        enrichment_dic = {}
        for line in filein_csv:
            enrichment_dic[line[0]]=float(line[1])   
except:
    sys.exit('Cannot find the following enrichment file:\n' + file_in)

###############################################################################
# Get phenotypes and save files
if record_dic['guide_type']=='single':
    gene2enrichment, gene2elementname, gene2medphe, gene2tscore, gene2bootstrap, gene2mwu = calculatePhe(enrichment_dic,record_dic,args)
elif record_dic['guide_type']=='double':
    front_e_singleko, back_e_singleko, front_e_singleko_medphe, back_e_singleko_medphe,\
    double_name_list, double_expphe_list, double_obsphe_list, rawGI_list, normGI_list, neg_name_list, neg_expphe_list, neg_obsphe_list,\
    o_genepair2expphe, o_genepair2obsphe, o_genepair2normgi,\
    o_genepair2medexpphe, o_genepair2medobsphe, o_genepair2mednormgi,\
    genepair2expphe, genepair2obsphe, genepair2normgi,\
    genepair2medexpphe, genepair2medobsphe, genepair2mednormgi,\
    o_gene2bootstrap, o_genepair2tscore, o_genepair2mwu, gene2bootstrap, genepair2tscore, genepair2mwu = calculatePhe(enrichment_dic,record_dic,args)
    
    # write normGI values for each sgRNAs
    file_in = args.file_head + '_normGI.csv'
    with open(file_in,'wb') as temp:
        tempcsv = csv.writer(temp,delimiter=',')
        for i,pgRNA in enumerate(double_name_list):
            tempcsv.writerow([pgRNA,normGI_list[i]])