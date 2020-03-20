#####################################################################
#####################################################################

## (5) Combine two replicates to calculate gene phenotypes

## python CombineReplicates.py <sgRNA enrichment csv file> <output file>

## ex) python CombineReplicates.py ../results/Batch_retest/T0_vs_Day21_23_rep1 ../results/Batch_retest/T0_vs_Day21_23_rep2 ../results/Batch_retest/T0_vs_Day21_23_combo

#####################################################################
#####################################################################


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
parser.add_argument('file_head1', help='file head for replicate1 enrichment and record input file', type=str)
parser.add_argument('file_head2', help='file head for replicate2 enrichment and record input file', type=str)
parser.add_argument('file_head', help='file head for output file', type=str)

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

file_in = args.file_head1 + '_record.txt'
try:
    with open(file_in, 'r') as in_open:
        print('Record file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter='\t')
        record_dic1 = {}
        for line in filein_csv:
            record_dic1[line[0]]=line[1]   
except:
    sys.exit('Cannot find the following record file:\n' + file_in)

file_in = args.file_head1 + '_enrichment.csv'
try:
    with open(file_in, 'r') as in_open:
        print('Enrichment file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter=',')
        enrichment_dic1 = {}
        for line in filein_csv:
            enrichment_dic1[line[0]]=float(line[1])   
except:
    sys.exit('Cannot find the following enrichment file:\n' + file_in)

file_in = args.file_head2 + '_record.txt'
try:
    with open(file_in, 'r') as in_open:
        print('Record file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter='\t')
        record_dic2 = {}
        for line in filein_csv:
            record_dic2[line[0]]=line[1]   
except:
    sys.exit('Cannot find the following record file:\n' + file_in)

file_in = args.file_head2 + '_enrichment.csv'
try:
    with open(file_in, 'r') as in_open:
        print('Enrichment file found: '+file_in+' found')
        filein_csv = csv.reader(in_open,delimiter=',')
        enrichment_dic2 = {}
        for line in filein_csv:
            enrichment_dic2[line[0]]=float(line[1])   
except:
    sys.exit('Cannot find the following enrichment file:\n' + file_in)    

if record_dic1['guide_type']=='single':
    
    file_in = args.file_head1 + '_phenotype.csv'
    try:
        with open(file_in, 'r') as in_open:
            print('Phenotype file found: '+file_in+' found')
            filein_csv = csv.reader(in_open,delimiter=',')
            
            phenotype_dic1 = {}
            tscore_dic1 = {}
            bootstrap_dic1 = {}
            mwu_dic1 = {}
            
            for i,line in enumerate(filein_csv):
                if i==0:
                    continue
                phenotype_dic1[line[2]]=float(line[8])
                try:
                    tscore_dic1[line[2]]=float(line[9])
                except:
                    pass
                try:
                    bootstrap_dic1[line[2]]=float(line[10])
                except:
                    pass
                try:
                    mwu_dic1[line[2]]=float(line[11])
                except:
                    pass
    except:
        sys.exit('Cannot find the following phenotype file:\n' + file_in)    
    
    file_in = args.file_head2 + '_phenotype.csv'
    try:
        with open(file_in, 'r') as in_open:
            print('Phenotype file found: '+file_in+' found')
            filein_csv = csv.reader(in_open,delimiter=',')
            
            phenotype_dic2 = {}
            tscore_dic2 = {}
            bootstrap_dic2 = {}
            mwu_dic2 = {}
            
            for i,line in enumerate(filein_csv):
                if i==0:
                    continue                
                phenotype_dic2[line[2]]=float(line[8])
                try:
                    tscore_dic2[line[2]]=float(line[9])
                except:
                    pass
                try:
                    bootstrap_dic2[line[2]]=float(line[10])
                except:
                    pass
                try:
                    mwu_dic2[line[2]]=float(line[11])
                except:
                    pass
    except:
        sys.exit('Cannot find the following phenotype file:\n' + file_in)    

elif record_dic1['guide_type']=='double':

    file_in = args.file_head1 + '_phenotype.csv'
    try:
        with open(file_in, 'r') as in_open:
            print('Phenotype file found: '+file_in+' found')
            filein_csv = csv.reader(in_open,delimiter=',')
            
            phenotype_dic1 = {}
            normgi_dic1 = {}
            tscore_dic1 = {}
            bootstrap_dic1 = {}
            mwu_dic1 = {}
            
            for i,line in enumerate(filein_csv):
                if i==0:
                    continue
                phenotype_dic1[line[0]]=float(line[10])
                normgi_dic1[line[0]]=float(line[12])
                try:                
                    tscore_dic1[line[0]]=float(line[13])
                except:
                    pass
                try:
                    bootstrap_dic1[line[0]]=float(line[14])
                except:
                    pass
                try:
                    mwu_dic1[line[0]]=float(line[15])
                except:
                    pass
    except:
        sys.exit('Cannot find the following phenotype file:\n' + file_in)    
    
    file_in = args.file_head2 + '_phenotype.csv'
    try:
        with open(file_in, 'r') as in_open:
            print('Phenotype file found: '+file_in+' found')
            filein_csv = csv.reader(in_open,delimiter=',')
            
            phenotype_dic2 = {}
            normgi_dic2 = {}
            tscore_dic2 = {}
            bootstrap_dic2 = {}
            mwu_dic2 = {}
            
            for i,line in enumerate(filein_csv):
                if i==0:
                    continue
                phenotype_dic2[line[0]]=float(line[10])
                normgi_dic2[line[0]]=float(line[12])
                try:               
                    tscore_dic2[line[0]]=float(line[13])
                except:
                    pass
                try:
                    bootstrap_dic2[line[0]]=float(line[14])
                except:
                    pass
                try:
                    mwu_dic2[line[0]]=float(line[15])
                except:
                    pass
    except:
        sys.exit('Cannot find the following phenotype file:\n' + file_in)
        
## Compare two record files to see if you have good replicate files

if record_dic1['guide_type']!=record_dic2['guide_type']:
    sys.exit('guide type is different between replicates')
    
if record_dic1['neg_ctrl']!=record_dic2['neg_ctrl']:
    sys.exit('neg_ctrl head is different between replicates')

if record_dic1['neg_idx']!=record_dic2['neg_idx']:
    sys.exit('neg_idx is different between replicates')

if record_dic1['split_double']!=record_dic2['split_double']:
    sys.exit('split_double is different between replicates')

if record_dic1['split_single']!=record_dic2['split_single']:
    sys.exit('split_single is different between replicates')

if record_dic1['gene_idx']!=record_dic2['gene_idx']:
    sys.exit('gene_idx is different between replicates')
    
if record_dic1['count_thresh']!=record_dic2['count_thresh']:
    sys.exit('count_thresh is different between replicates')    

if record_dic1['replace_zero']!=record_dic2['replace_zero']:
    sys.exit('replace_zero is different between replicates')        
    
###############################################################################
# Combine enrichment file to make new enrichment files
if record_dic1['guide_type']=='single':
    enrichment_dic = {}
    
    for key, value in enrichment_dic1.iteritems():
        newkey = key+record_dic1['split_single']+'rep1'
        enrichment_dic[newkey]=value
        
    for key, value in enrichment_dic2.iteritems():
        newkey = key+record_dic2['split_single']+'rep2'
        enrichment_dic[newkey]=value
    
elif record_dic1['guide_type']=='double':
    enrichment_dic = {}

    for key, value in enrichment_dic1.iteritems():
        front_e = key.split(record_dic1['split_double'])[0]
        back_e = key.split(record_dic1['split_double'])[1]
        
        new_front_e = front_e + record_dic1['split_single']+'rep1'
        new_back_e = back_e + record_dic1['split_single']+'rep1'
        
        newkey = new_front_e + record_dic1['split_double'] + new_back_e
        enrichment_dic[newkey]=value
        
    for key, value in enrichment_dic2.iteritems():
        front_e = key.split(record_dic2['split_double'])[0]
        back_e = key.split(record_dic2['split_double'])[1]
        
        new_front_e = front_e + record_dic2['split_single']+'rep2'
        new_back_e = back_e + record_dic2['split_single']+'rep2'
        
        newkey = new_front_e + record_dic2['split_double'] + new_back_e
        enrichment_dic[newkey]=value

###############################################################################
# Get phenotypes and save files
if record_dic1['guide_type']=='single':
    gene2enrichment, gene2elementname, gene2medphe, gene2tscore, gene2bootstrap, gene2mwu = calculatePhe(enrichment_dic,record_dic1,args)
elif record_dic1['guide_type']=='double':
    front_e_singleko, back_e_singleko, front_e_singleko_medphe, back_e_singleko_medphe,\
    double_name_list, double_expphe_list, double_obsphe_list, rawGI_list, normGI_list, neg_name_list, neg_expphe_list, neg_obsphe_list,\
    o_genepair2expphe, o_genepair2obsphe, o_genepair2normgi,\
    o_genepair2medexpphe, o_genepair2medobsphe, o_genepair2mednormgi,\
    genepair2expphe, genepair2obsphe, genepair2normgi,\
    genepair2medexpphe, genepair2medobsphe, genepair2mednormgi,\
    o_gene2bootstrap, o_genepair2tscore, o_genepair2mwu, gene2bootstrap, genepair2tscore, genepair2mwu = calculatePhe(enrichment_dic,record_dic1,args)
    
    # write normGI values for each sgRNAs
    file_in = args.file_head + '_rep1vs2_pgRNA_normGI.csv'
    with open(file_in,'wb') as temp:
        tempcsv = csv.writer(temp,delimiter=',')
        for i,pgRNA in enumerate(double_name_list):
            tempcsv.writerow([pgRNA,normGI_list[i]])     
###############################################################################
# Compare replicates

if record_dic1['guide_type']=='single':
    
    # phenotype repl vs rep2
    xx=[]; yy=[];
    for key in phenotype_dic1:
        if key in phenotype_dic2:
            xx.append(phenotype_dic1[key])
            yy.append(phenotype_dic2[key])

    # remove nan
    aa = np.isnan(xx)
    bb = np.isnan(yy)
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]
    
    # remove inf
    aa = (np.array(np.abs(xx))==float("inf"))
    bb = (np.array(np.abs(yy))==float("inf"))
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]

    Max = max(xx+yy)
    Min = min(xx+yy)
    Slope = np.polyfit(xx,yy,1)[0]
    Intercept = np.polyfit(xx,yy,1)[1]
                    
    plt.figure(figsize=(10,10))
    plt.plot(xx,yy,'k.',markersize=5,alpha=0.5,rasterized=True)
    Buff= (Max-Min)/5.0
    plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
    plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
    plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
    plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
    
    PearsonCoeff=st.stats.pearsonr(xx,yy)
    plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
    
    plt.xlim([Min-Buff,Max+Buff])
    plt.ylim([Min-Buff,Max+Buff])
    plt.xlabel('phenotype rep1 (pZ)')
    plt.ylabel('phenotype rep2 (pZ)')
    plt.tight_layout()
    plt.legend(loc='best')        
    
    file_in = args.file_head + '_rep1vs2_phenotype.pdf'        
    plt.savefig(file_in,format='PDF',transparent=True,dpi=400)        
    plt.close()
    
    # tscore repl vs rep2        
    xx=[]; yy=[];
    for key in tscore_dic1:
        if key in tscore_dic2:
            xx.append(tscore_dic1[key])
            yy.append(tscore_dic2[key])

    # remove nan
    aa = np.isnan(xx)
    bb = np.isnan(yy)
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]
    
    # remove inf
    aa = (np.array(np.abs(xx))==float("inf"))
    bb = (np.array(np.abs(yy))==float("inf"))
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]
                        
    if len(xx)!=0:
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]        
                                
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'k.',markersize=5,alpha=0.5,rasterized=True)
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('tscore rep1')
        plt.ylabel('tscore rep2')
        plt.tight_layout()
        plt.legend(loc='best')                
        file_in = args.file_head + '_rep1vs2_tscore.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)        
        plt.close()

    # bootstrap pval repl vs rep2        
    xx=[]; yy=[];
    for key in bootstrap_dic1:
        if key in bootstrap_dic2:

            if bootstrap_dic1[key]==0:
                xv=-1*np.log10(0.1/args.bootstrap_repeat)*np.sign(phenotype_dic1[key])
            else:
                xv=-1*np.log10(bootstrap_dic1[key])*np.sign(phenotype_dic1[key])
                
            if bootstrap_dic2[key]==0:
                yv=-1*np.log10(0.1/args.bootstrap_repeat)*np.sign(phenotype_dic2[key])
            else:
                yv=-1*np.log10(bootstrap_dic2[key])*np.sign(phenotype_dic2[key])
            
            xx.append(xv)
            yy.append(yv)

    # remove nan
    aa = np.isnan(xx)
    bb = np.isnan(yy)
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]
    
    # remove inf
    aa = (np.array(np.abs(xx))==float("inf"))
    bb = (np.array(np.abs(yy))==float("inf"))
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]

    if len(xx)!=0:    
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
                        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'k.',markersize=5,alpha=0.5,rasterized=True)
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('signed log10 pval rep1')
        plt.ylabel('signed log10 pval rep2')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = args.file_head + '_rep1vs2_bootstrap_signed_log10pval.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()

    # mwu pval repl vs rep2        
    xx=[]; yy=[];
    for key in mwu_dic1:
        if key in mwu_dic2:

            xv=-1*np.log10(mwu_dic1[key])*np.sign(phenotype_dic1[key])
            yv=-1*np.log10(mwu_dic2[key])*np.sign(phenotype_dic2[key])
            
    # remove nan
    aa = np.isnan(xx)
    bb = np.isnan(yy)
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]
    
    # remove inf
    aa = (np.array(np.abs(xx))==float("inf"))
    bb = (np.array(np.abs(yy))==float("inf"))
    cc = np.logical_or(aa,bb)
    
    xx = np.array(xx)[~cc]
    yy = np.array(yy)[~cc]

    if len(xx)!=0:    
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
                        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'k.',markersize=5,alpha=0.5,rasterized=True)
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('signed log10 pval rep1')
        plt.ylabel('signed log10 pval rep2')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = args.file_head + '_rep1vs2_mwu_signed_log10pval.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()               
                                
if record_dic1['guide_type']=='double':
   
    # phenotype rep1 vs rep2
    xx=[]; yy=[]; 
    for key in phenotype_dic1:
        if key in phenotype_dic2:
            xx.append(phenotype_dic1[key])
            yy.append(phenotype_dic2[key])

    Max = max(xx+yy)
    Min = min(xx+yy)
    Slope = np.polyfit(xx,yy,1)[0]
    Intercept = np.polyfit(xx,yy,1)[1]
                    
    plt.figure(figsize=(10,10))
    plt.plot(xx,yy,'k.',markersize=5,alpha=0.5,rasterized=True)
    Buff= (Max-Min)/5.0
    plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
    plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
    plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
    plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
    
    PearsonCoeff=st.stats.pearsonr(xx,yy)
    plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
    
    plt.xlim([Min-Buff,Max+Buff])
    plt.ylim([Min-Buff,Max+Buff])
    plt.xlabel('phenotype rep1 (pZ)')
    plt.ylabel('phenotype rep2 (pZ)')
    plt.tight_layout()
    plt.legend(loc='best')        
    
    file_in = args.file_head + '_rep1vs2_phenotype.pdf'        
    plt.savefig(file_in,format='PDF',transparent=True,dpi=400)        
    plt.close()

    # normgi rep1 vs rep2
    xx=[]; yy=[]; xn=[]; yn=[]; xnn=[]; ynn=[]; xs=[]; ys=[];
    for key in normgi_dic1:
        if key in normgi_dic2:
            xx.append(normgi_dic1[key])
            yy.append(normgi_dic2[key])
            
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==1:
                xn.append(normgi_dic1[key])
                yn.append(normgi_dic2[key])
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==2:
                xnn.append(normgi_dic1[key])
                ynn.append(normgi_dic2[key])
            if front_g==back_g:
                xs.append(normgi_dic1[key])
                ys.append(normgi_dic2[key])

    Max = max(xx+yy)
    Min = min(xx+yy)
    Slope = np.polyfit(xx,yy,1)[0]
    Intercept = np.polyfit(xx,yy,1)[1]
                
    plt.figure(figsize=(10,10))
    plt.plot(xx,yy,'ko',markersize=5,alpha=0.5,rasterized=True,label='all')
    plt.plot(xn,yn,'yo',markersize=10,alpha=1,label='x_negctrl')
    plt.plot(xs,ys,'ro',markersize=5,alpha=1,label='x_x')    
    plt.plot(xnn,ynn,'co',markersize=10,alpha=1,label='negctrl_negctrl')    
    
    Buff= (Max-Min)/5.0
    plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
    plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
    plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
    plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
    
    PearsonCoeff=st.stats.pearsonr(xx,yy)
    plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
    
    plt.xlim([Min-Buff,Max+Buff])
    plt.ylim([Min-Buff,Max+Buff])
    plt.xlabel('normGI rep1')
    plt.ylabel('normGI rep2')
    plt.tight_layout()
    plt.legend(loc='best')        
    
    file_in = args.file_head + '_rep1vs2_normgi.pdf'        
    plt.savefig(file_in,format='PDF',transparent=True,dpi=400)        
    plt.close()

    # tscore rep1 vs rep2                                                                                                
    xx=[]; yy=[]; xn=[]; yn=[]; xnn=[]; ynn=[]; xs=[]; ys=[];
    for key in tscore_dic1:
        if key in tscore_dic2:
            xx.append(tscore_dic1[key])
            yy.append(tscore_dic2[key])
            
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==1:
                xn.append(tscore_dic1[key])
                yn.append(tscore_dic2[key])
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==2:
                xnn.append(tscore_dic1[key])
                ynn.append(tscore_dic2[key])
            if front_g==back_g:
                xs.append(tscore_dic1[key])
                ys.append(tscore_dic2[key])                         
            
    if len(xx)!=0:
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]        
                                
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'ko',markersize=5,alpha=0.5,rasterized=True,label='all')
        plt.plot(xn,yn,'yo',markersize=10,alpha=1,label='x_negctrl')
        plt.plot(xs,ys,'ro',markersize=5,alpha=1,label='x_x')    
        plt.plot(xnn,ynn,'co',markersize=10,alpha=1,label='negctrl_negctrl')
        
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('tscore rep1')
        plt.ylabel('tscore rep2')
        plt.tight_layout()
        plt.legend(loc='best')                
        file_in = args.file_head + '_rep1vs2_tscore.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)        
        plt.close()

    # bootstrap rep1 vs rep2
    xx=[]; yy=[]; xn=[]; yn=[]; xnn=[]; ynn=[]; xs=[]; ys=[];
    for key in bootstrap_dic1:
        if key in bootstrap_dic2:

            if bootstrap_dic1[key]==0:
                xv=-1*np.log10(0.1/args.bootstrap_repeat)*np.sign(normgi_dic1[key])
            else:
                xv=-1*np.log10(bootstrap_dic1[key])*np.sign(normgi_dic1[key])
                
            if bootstrap_dic2[key]==0:
                yv=-1*np.log10(0.1/args.bootstrap_repeat)*np.sign(normgi_dic2[key])
            else:
                yv=-1*np.log10(bootstrap_dic2[key])*np.sign(normgi_dic2[key])
            
            xx.append(xv)
            yy.append(yv)
            
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==1:
                xn.append(xv)
                yn.append(yv)
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==2:
                xnn.append(xv)
                ynn.append(yv)   
            if front_g==back_g:
                xs.append(xv)
                ys.append(yv)                                            

    if len(xx)!=0:    
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
                        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'ko',markersize=5,alpha=0.5,rasterized=True,label='all')
        plt.plot(xn,yn,'yo',markersize=10,alpha=1,label='x_negctrl')
        plt.plot(xs,ys,'ro',markersize=5,alpha=1,label='x_x')    
        plt.plot(xnn,ynn,'co',markersize=10,alpha=1,label='negctrl_negctrl')    
        
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('signed log10 pval rep1')
        plt.ylabel('signed log10 pval rep2')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = args.file_head + '_rep1vs2_bootsrap_signed_log10pval.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()

    # mwu rep1 vs rep2
    xx=[]; yy=[]; xn=[]; yn=[]; xnn=[]; ynn=[]; xs=[]; ys=[];
    for key in mwu_dic1:
        if key in mwu_dic2:
            xv=-1*np.log10(mwu_dic1[key])*np.sign(normgi_dic1[key])
            yv=-1*np.log10(mwu_dic2[key])*np.sign(normgi_dic2[key])
            
            xx.append(xv)
            yy.append(yv)
            
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==1:
                xn.append(xv)
                yn.append(yv)
            if (front_g==record_dic1['neg_ctrl'])+(back_g==record_dic1['neg_ctrl'])==2:
                xnn.append(xv)
                ynn.append(yv)   
            if front_g==back_g:
                xs.append(xv)
                ys.append(yv)                                            

    if len(xx)!=0:    
        Max = max(xx+yy)
        Min = min(xx+yy)
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
                        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'ko',markersize=5,alpha=0.5,rasterized=True,label='all')
        plt.plot(xn,yn,'yo',markersize=10,alpha=1,label='x_negctrl')
        plt.plot(xs,ys,'ro',markersize=5,alpha=1,label='x_x')    
        plt.plot(xnn,ynn,'co',markersize=10,alpha=1,label='negctrl_negctrl')    
        
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('signed log10 pval rep1')
        plt.ylabel('signed log10 pval rep2')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = args.file_head + '_rep1vs2_mwu_signed_log10pval.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()
              