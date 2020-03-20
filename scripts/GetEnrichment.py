#####################################################################################
# Kyuho Han (Some codes are adapted from David Morgen's scripts)
# 10/10/2017
# Compare two count files to calculate normalized log fold enrichment of sgRNAs
#####################################################################################
# run GetEnrichment.py ../counts/20171009_SL_Screens/20171009_SL_Plas_pgDouble_Final.counts ../counts/20171009_SL_Screens/20171009_SL_D15_Unt1_pgDouble_Final.counts ../results/20171009_SL_Screens/Plas_vs_D15_Unt1
# run GetEnrichment.py ../counts/20171009_SL_Screens/20171009_SL_Plas_pgDouble_Final.counts ../counts/20171009_SL_Screens/20171009_SL_D15_Unt2_pgDouble_Final.counts ../results/20171009_SL_Screens/Plas_vs_D15_Unt2

# run GetEnrichment.py ../counts/20171031_KRAS_3D/T0_counts.csv ../counts/20171031_KRAS_3D/D19_2D_U1_counts.csv ../results/20171031_KRAS_3D/T0_vs_D19_2D_U1 -gt single -ct 10
# run GetEnrichment.py ../counts/20171208_Youngtae/D16_CONT_A_counts.csv ../counts/20171208_Youngtae/D16_IRA_A_counts.csv ../results/20171208_Youngtae/D16_CONT_A_vs_D16_IRA_A -gt single -ct 10

from __future__ import division

import matplotlib as mpl 
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
plt.close("all")

import sys
import csv
import os
import argparse
from FunctionLib import *

import numpy as np
import scipy.stats as st

#####################################################################################
# Version number
current_version = '0.1'

#####################################################################################
# Parse input
parser = argparse.ArgumentParser(description='Obtain normalized log fold enrichment of sgRNAs/shRNAs (pZ)')

# Non-optional arguments
parser.add_argument('unt_file', help='File for untreated counts', type=str)
parser.add_argument('trt_file', help='File for treated counts', type=str)
parser.add_argument('output', help='Name for output files', type=str)

# Options for element IDs
parser.add_argument('-gt', '--guide_type',help='Types of sgRNAs/shRNAs used. Default is double', type=str, default='double', choices=['single','double'])
parser.add_argument('-nc', '--neg_ctrl',help='Symbol used to denote negative controls. Default is safe.',type=str, default='safe')
parser.add_argument('-ni', '--neg_idx',help='Position of negname in parsed element name. Default is 1.',type=int, default=1)
parser.add_argument('-se', '--split_double',help='Delimiter separting two elements in double. Default is ____',type=str, default='____')
parser.add_argument('-sn', '--split_single',help='Delimiter for parsing single element name. Default is _',type=str, default='_')
parser.add_argument('-gi', '--gene_idx',help='Position of genename in parsed element name. Default is 1', type=int, default=1)

# Options for analysis
parser.add_argument('-ct', '--count_thresh', help='Read cutoff for small count numbers. Default is 50.',type=int, default=50)
parser.add_argument('-rz', '--replace_zero', help='Replace zero count with a certain number. Default is 1',type=int,default=1)
parser.add_argument('-sc', '--sign_correction', help='Sign correction for GI. Default is y', type=str, default='y', choices=['y','n'])

# Options for computation
parser.add_argument('-p', '--proccessors', dest='nums',help='Number of proccessors to use; default is 20', type=int,default=20)

args = parser.parse_args()

#####################################################################################
# Processes and checks input arguments

# Creates output file name
file_out = os.path.join(args.output)

# Checks if it can write to output
try:
    with open(file_out + '_record.txt', 'w') as out_open:
        pass
    os.remove(file_out + '_record.txt')
except:
    sys.exit('Cannot write to output file:\n'+file_out)

#####################################################################################
# Import count files, filter out sgRNAs with low threshold count, and calculate enrichment

print('Calculating enrichment of sgRNAs/shRNAs')
unt_count, trt_count, enrichment, guides_removed = CalEnrichment(args.unt_file,args.trt_file,args.count_thresh,args.replace_zero)

print('Total untreated count: '+str(sum(unt_count.values())))
print('Total treated count: '+str(sum(trt_count.values())))
print('Total elements detected: '+str(len(list(set(unt_count.keys()+trt_count.keys())))))
print('Total elements removed (under '+str(args.count_thresh)+' reads): '+str(len(guides_removed)))

#####################################################################################
# Plotting count files, unt vs trt 

plt.figure(figsize=(10,10))

xval = []; yval = []; xneg = []; yneg = [];
for key in unt_count:
    if key in trt_count:
        xval.append(np.log10(unt_count[key]))
        yval.append(np.log10(trt_count[key]))
        
        if args.guide_type == 'double':
            front_e = key.split(args.split_double)[0]
            back_e = key.split(args.split_double)[1]
            
            if front_e.split(args.split_single)[args.neg_idx]==args.neg_ctrl and back_e.split(args.split_single)[args.neg_idx]==args.neg_ctrl:
                xneg.append(np.log10(unt_count[key]))
                yneg.append(np.log10(trt_count[key]))
        elif args.guide_type == 'single':
            if key.split(args.split_single)[args.neg_idx]==args.neg_ctrl:
                xneg.append(np.log10(unt_count[key]))
                yneg.append(np.log10(trt_count[key]))

plt.plot(xval,yval,'k.',rasterized=True,label='all elements',alpha=0.2)
plt.plot(xneg,yneg,'r.',rasterized=True,label='neg ctrls',alpha=0.5)
plt.xlabel('unt count, log10')
plt.ylabel('trt count, log10')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(file_out+'_unt_vs_trt_count.pdf',format='PDF',dpi=400,transparent=True)
plt.close()

#####################################################################################
# Normalize enrichment by negative controls
norm_enrichment, neg_enrichment = NormEnrichment(enrichment,args.split_single,args.split_double,args.neg_ctrl,args.neg_idx,args.guide_type)

print('Negative control detected: '+ str(len(neg_enrichment)))

#####################################################################################
# Compare two different orientations for double
if args.guide_type=='double':
    xval=[]; yval=[]; labels=[];
    
    check_overlap = []
    for key in norm_enrichment:
        front_e=key.split(args.split_double)[0]
        back_e=key.split(args.split_double)[1]
        
        flipped_e = back_e + args.split_double + front_e
        check_overlap.append(key)
        
        if (flipped_e in norm_enrichment) and (flipped_e not in check_overlap):
            xval.append(norm_enrichment[key])
            yval.append(norm_enrichment[flipped_e])
            labels.append(key)
    
    Max = max(xval+yval)
    Min = min(xval+yval)    
    
    Slope = np.polyfit(xval,yval,1)[0]
    Intercept = np.polyfit(xval,yval,1)[1]

    plt.figure(figsize=(10,10))
    plt.plot(xval,yval,'k.',markersize=5,alpha=.2,rasterized=True)
    Buff= (Max-Min)/5.0
    plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
    plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
    plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
    plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])

    PearsonCoeff=st.stats.pearsonr(xval,yval)
    plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
    plt.legend(loc='best')
    
    plt.xlim([Min-Buff,Max+Buff])
    plt.ylim([Min-Buff,Max+Buff])
    plt.xlabel('AB_double_element (pZ)')
    plt.ylabel('BA_double_element (pZ)')
    plt.tight_layout()
    plt.savefig(file_out+'_ABBA_double_sgRNA_phe_plot.pdf',format='PDF',transparent=True,dpi=400)
    plt.close()     
    
#####################################################################################
# Write output files : enrichment values of sgRNA/shRNAs, record file
with open(file_out + '_enrichment.csv', 'wb') as out_open:
    out_csv = csv.writer(out_open, delimiter=',')
    for key, value in norm_enrichment.iteritems():
        out_csv.writerow([key,value])

with open(file_out + '_record.txt', 'wb') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['GetEnrichment.py', current_version])
    rec_csv.writerow(['unt_file', args.unt_file])
    rec_csv.writerow(['trt_file', args.trt_file])
    rec_csv.writerow(['guide_type', args.guide_type])
    rec_csv.writerow(['neg_ctrl', args.neg_ctrl])
    rec_csv.writerow(['neg_idx', args.neg_idx])
    rec_csv.writerow(['split_double', args.split_double])
    rec_csv.writerow(['split_single', args.split_single])
    rec_csv.writerow(['gene_idx', args.gene_idx])
    rec_csv.writerow(['count_thresh', args.count_thresh])
    rec_csv.writerow(['replace_zero', args.replace_zero]) 
    rec_csv.writerow(['sign_correction', args.sign_correction]) 
        
with open('../records/' + args.output.split('/')[-1] + '_record.txt', 'wb') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['GetEnrichment.py', current_version])
    rec_csv.writerow(['unt_file', args.unt_file])
    rec_csv.writerow(['trt_file', args.trt_file])
    rec_csv.writerow(['guide_type', args.guide_type])
    rec_csv.writerow(['neg_ctrl', args.neg_ctrl])
    rec_csv.writerow(['neg_idx', args.neg_idx])
    rec_csv.writerow(['split_double', args.split_double])
    rec_csv.writerow(['split_single', args.split_single])
    rec_csv.writerow(['gene_idx', args.gene_idx])
    rec_csv.writerow(['count_thresh', args.count_thresh])
    rec_csv.writerow(['replace_zero', args.replace_zero])     
    rec_csv.writerow(['sign_correction', args.sign_correction]) 