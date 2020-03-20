###############################################################################
# Kyuho Han # Some codes are adapted from David Morgen's scripts 10/10/2017
###############################################################################

from __future__ import division

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
plt.close("all")
import matplotlib.cm as cm

import numpy as np
import csv
from collections import defaultdict
import os
import scipy.misc
import scipy.stats as st
import scipy.stats.mstats as ms

from statsmodels.sandbox.stats.multicomp import multipletests

import sys
import random
import re
import timeit
from scipy.interpolate import interp1d
from patsy import dmatrix
# import pandas as pd

###############################################################################

def CalEnrichment(unt_file, trt_file, count_thresh, replace_zero):
    unt_count = {}
    trt_count = {}
    
    with open(unt_file, 'rU') as unt_open:
        dialect = csv.Sniffer().sniff(unt_open.read(1024), delimiters='\t ,')
        unt_open.seek(0)
        unt_csv = csv.reader(unt_open, dialect)
        
        for line in unt_csv:
            unt_count[line[0]]=int(line[1])
            
    with open(trt_file, 'rU') as trt_open:
        dialect = csv.Sniffer().sniff(trt_open.read(1024), delimiters='\t ,')
        trt_open.seek(0)
        trt_csv = csv.reader(trt_open, dialect)
        
        for line in trt_csv:
            trt_count[line[0]]=int(line[1])
            
    unique_elements = list(set(unt_count.keys() + trt_count.keys()))
    
    enrichment = {}
    guides_removed = {}
    
    totaluntcount = np.sum(unt_count.values())
    totaltrtcount = np.sum(trt_count.values())
    countratio = totaluntcount / totaltrtcount
    
    for element in unique_elements:        
        if element in unt_count:
            uc = unt_count[element]
        else:
            uc = replace_zero
            
        if uc == 0:
            uc = replace_zero

        if element in trt_count:
            tc = trt_count[element]
        else:
            tc = replace_zero
            
        if tc == 0:
            tc = replace_zero
            
        if uc >= count_thresh: # or tc >= count_thresh:
            enrichment[element] = np.log2(countratio*tc/uc)
        else:
            guides_removed[element]=[uc, tc]
            
    return unt_count, trt_count, enrichment, guides_removed
            
def NormEnrichment(enrichment, split_single, split_double, neg_ctrl, neg_idx, guide_type):
    if guide_type=='double':
        neg_enrichment = {}
        for key, value in enrichment.iteritems():
            front_element = key.split(split_double)[0]
            back_element = key.split(split_double)[1]
            
            if front_element.split(split_single)[neg_idx] == neg_ctrl and back_element.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue
        
        neglog2 = np.median(neg_enrichment.values())
        
        # subtract neglog2
        norm1_enrichment = {}    
        for key, value in enrichment.iteritems():
            norm1_enrichment[key] = value - neglog2
    
        neg_enrichment = {}
        for key, value in norm1_enrichment.iteritems():
            front_element = key.split(split_double)[0]
            back_element = key.split(split_double)[1]
            
            if front_element.split(split_single)[neg_idx] == neg_ctrl and back_element.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue
        
        # divide enrichment by std of negctrl
        norm2_enrichment = {}
        sigma = np.std(neg_enrichment.values())
        for key, value in norm1_enrichment.iteritems():
            norm2_enrichment[key]=value/sigma

        neg_enrichment = {}
        for key, value in norm2_enrichment.iteritems():
            front_element = key.split(split_double)[0]
            back_element = key.split(split_double)[1]
            
            if front_element.split(split_single)[neg_idx] == neg_ctrl and back_element.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue            
    
    elif guide_type=='single':
        neg_enrichment = {}
        for key, value in enrichment.iteritems():        
            if key.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue
        
        neglog2 = np.median(neg_enrichment.values())
        
        # subtract neglog2
        norm1_enrichment = {}    
        for key, value in enrichment.iteritems():
            norm1_enrichment[key] = value - neglog2
    
        neg_enrichment = {}
        for key, value in norm1_enrichment.iteritems():
            if key.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue
        
        # divide enrichment by std of negctrl
        norm2_enrichment = {}
        sigma = np.std(neg_enrichment.values())
        for key, value in norm1_enrichment.iteritems():
            norm2_enrichment[key]=value/sigma 

        neg_enrichment = {}
        for key, value in norm2_enrichment.iteritems():
            if key.split(split_single)[neg_idx] == neg_ctrl:
                neg_enrichment[key]=value
            else:
                continue            
    
    return norm2_enrichment, neg_enrichment
            
###############################################################################
# Function which retrieves info about genes (David Morgens' code)

def retrieveInfo(ref_base='../LUT', mouse=False):

    '''
    Retrieves gene info for the screen type. Location of reference
    files can be changed, defaults to nearby GenRef folder.
    '''

    # Finds info files downloaded from NCBI
    org_file_human = os.path.join(ref_base, 'Homo_sapiens.gene_info')
    org_file_mouse = os.path.join(ref_base, 'Mus_musculus.gene_info')

    # Custom Ensemble ID to gene name file
    ens_file = os.path.join(ref_base, 'ensRef.csv')

    geneID2Name = defaultdict(lambda: 'N/A')
    geneID2Info = defaultdict(lambda: 'N/A')
    geneName2ID = defaultdict(lambda: 'N/A')
    geneEns2Name = defaultdict(lambda: 'N/A')
    geneName2Ens = defaultdict(lambda: 'N/A')

    # Reads in Ensemble data
    try:
        with open(ens_file, 'r') as ens_open:

            ens_csv = csv.reader(ens_open, delimiter=',')

            for line in ens_csv:
                geneEns2Name[line[1]] = line[0]
                geneName2Ens[line[0]] = line[1]

    except IOError:
        print('Ensembl information file not found.\n'
                + 'Use -r to change file location')

    # Reads in Mouse data
    try:
        with open(org_file_mouse, 'r') as org_open:

            org_csv = csv.reader(org_open, delimiter='\t')
            org_csv.next()  # Skips header

            for line in org_csv:
                # Entrez
                geneID2Name[line[1]] = line[2]
                geneID2Info[line[1]] = line[8]
                geneName2ID[line[2]] = line[1]

    except IOError:
        print('Mouse information file not found.\n'
                + 'Use -r to change file location')

    if not mouse:
        # Reads in Human data
        try:
            with open(org_file_human, 'r') as org_open:

                org_csv = csv.reader(org_open, delimiter='\t')
                org_csv.next()  # Skips header

                # For each line in file, save that gene information
                for line in org_csv:

                    geneID2Name[line[1]] = line[2]
                    geneID2Info[line[1]] = line[8]
                    geneName2ID[line[2]] = line[1]


        except IOError:
            print('Human information file not found.\n'
                    + 'Use -r to change file location')

    return geneID2Name, geneID2Info, geneName2ID, geneEns2Name, geneName2Ens

###############################################################################
# Retreives GO information (David Morgens' code)

def retrieveGO(ref_base='../LUT'):
    '''
    Returns GO component, process, and function data by geneID.
    '''

    go_file = os.path.join(ref_base, 'gene2go')

    # Stores as dictionary of strings. 
    geneID2Comp = defaultdict(str)
    geneID2Proc = defaultdict(str)
    geneID2Fun = defaultdict(str)

    # Checks that file exists, if not returns empty dictionaries
    if not os.path.isfile(go_file):

        print('GO reference file not found; use -r to change file location')
        return geneID2Comp, geneID2Proc, geneID2Fun

    # Reads in GO data
    with open(go_file, 'r') as go_open:
        go_csv = csv.reader(go_open, delimiter='\t')

        for line in go_csv:

            # Checks that line is correct length
            if len(line) == 8:

                # Skips NOT GO terms
                if line[4] == 'NOT':
                    continue

                if line[7] == 'Component':
                    geneID2Comp[line[1]] += line[5] + '|'

                elif line[7] == 'Process':
                    geneID2Proc[line[1]] += line[5] + '|'

                elif line[7] == 'Function':
                    geneID2Fun[line[1]] += line[5] + '|'

    geneID2Comp.default_factory = lambda: 'None'
    geneID2Proc.default_factory = lambda: 'None'
    geneID2Fun.default_factory = lambda: 'None'

    return geneID2Comp, geneID2Proc, geneID2Fun
    
###############################################################################
# CalculatePhe

#def sample(data):
#    sample = [random.choice(data) for _ in xrange(len(data))]
#    return sample                
                
#def bootstrap_test(value,func=np.mean,repeat=10000):    
#    output = []
#    
#    for i in xrange(repeat):
#        output.append(func(np.random.choice(value,size=len(value))))
#    return np.array(output)
#
def bootstrap_test(value,func=np.mean,repeat=10000):    
    output = np.zeros((repeat,1))
    
    for i in xrange(repeat):
        output[i,0]=func(np.random.choice(value,size=len(value)))
    return output    
    
    
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N     

# def CalculateQvalFromPval(PvalDict,negId1,negId2,outbase,FDR=0.05):
# #Based on John D.Storey, Robert Tibashirani, PNAS, 2003 (statistical significance for genomewide studies)
# #FDR = PEI * m * t / S(t)
# # PEI: proportion of trully null (m0 / m) (unknown...)
# # m : total number of P values (ex: total number of samples, 20,000 genes)
# # t : p-val cutoff (ex: 0.05)
# # S(t) : the number of significant features when cutoff p-val is 't'
# 
#     #PvalDict = BetaGeneRawGIMWULMA_rep1n2
#     #negId1 = NegId1
#     #negId2 = NegId2
#     #outbase = OutBase
#     #FDR = 0.05    
# 
#     Pval= sorted(PvalDict.values())
#     bin_num = float(20)    
#     Lambda = Pval/max(Pval)
#     hist,edge = np.histogram(Lambda,bins=bin_num,range=(0,1))
#     num_norm = np.sum(hist)/float(np.sum(hist!=0))                                                                       
#     
#     #1) P value distribution
#     plt.figure(figsize=(18,12))
#     
#     plt.subplot(2,3,1)
#     plt.bar(edge[:-1], hist/num_norm, width=0.05)
#     plt.plot([0,1],[1,1],'k-.',alpha=0.5,label='expected density\nif all genes were null')                                              
#     plt.title('Density hitgoram of p values')
#     plt.xlabel('Lambda (Pval rescaled between 0 and 1)')
#     plt.ylabel('Density')
#     plt.xlim((0,1))
#     plt.legend(loc='best',fontsize=12)
#     
#     #2) Get PEI ( PEI : proportion of trully null : m0/m)
#     hist,edge = np.histogram(Lambda,bins=bin_num*2,range=(0,1))
#     num_norm = np.sum(hist)/float(np.sum(hist!=0))                                                                       
#     
#     XX = edge[1:]
#     YY = hist/num_norm
#     
#     YYY = dmatrix("cr(x, df=3) -1", {"x": XX})
#     
#     IDX1=int(len(XX)/6)
#     IDX2= IDX1*3
#     IDX3= IDX1*5
#     
#     B=np.array([list(YY)[IDX1],list(YY)[IDX2],list(YY)[IDX3]])
#     
#     plt.subplot(2,3,2)
#     
#     plt.plot(XX[IDX1:],YY[IDX1:],'ko',markersize=7)                                                                  
#     
#     plt.plot(XX,np.dot(YYY,B),'k-',alpha=0.5,label='natural cubic spine fit')
#     plt.xlim((-0.1,1.1))
#     plt.ylim((0,max(YY[IDX1:])*1.1))
#     plt.ylabel('PEI (proportion of truly null)')
#     plt.xlabel('Lambda (Pval rescaled from 0 to 1)')                                                                                                                                                                                                      
#     plt.title('PEI versus Lambda')
#     plt.legend(loc='best',fontsize=12)
#     plt.text(0.3,np.dot(YYY,B)[-1]*0.7,'PEI at (lambda=1) ='+str(np.dot(YYY,B)[-1])[0:4])
#  
#      #3) P-val Versus Q-val    
#     PEI = np.dot(YYY,B)[-1]
#     M = len(PvalDict)
#     QvalDict={}
#     
#     Pval_List=[]
#     Qval_List=[]
#     Key_List=[]
#     TotalNegNum=0
#     
#     for key in PvalDict:        
#         Key_List.append(key)
#         if negId1 in key or negId2 in key:
#             TotalNegNum+=1
#             
#         T = PvalDict[key]
#         Pval_List.append(T)
#         
#         S = np.sum(Pval<=T)
#     
#         Q = PEI * M * T / S
#         QvalDict[key] = Q
#         Qval_List.append(Q)
#     
#     IDX = np.argsort(Pval_List)
#     KeySel_List = np.array(Key_List)[np.array(Qval_List)<FDR]
#     SelNegNum=0
#     
#     for key in KeySel_List:
#         if negId1 in key or negId2 in key:
#             SelNegNum+=1
#         
#     Min_IDX = np.argmin(abs(np.array(Qval_List) - FDR))                                            
#     PVALatFDR = np.array(Pval_List)[Min_IDX]                                                                                            
#                                                                                                                                                                                                 
#     plt.subplot(2,3,3)
#     plt.plot(np.array(Pval_List)[IDX],np.array(Qval_List)[IDX],'r.',markersize=2,rasterized=True)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#     plt.xlabel('P values')
#     plt.ylabel('Q values')
#     plt.plot([0,PVALatFDR],[FDR,FDR],'k-',alpha=0.5)
#     plt.plot([PVALatFDR,PVALatFDR],[0,FDR],'k-',alpha=0.5)
#     plt.text(max(0.1, min(PVALatFDR*3,0.7)), max(0.1,min(FDR,0.9)), 'FDR: '+str(FDR)[0:4]+'\n'+'Pval: '+str(PVALatFDR)+'\nNo.of significants: '+str(len(KeySel_List)), fontsize=12) 
#     plt.title('P val VS Q val')                                       
#     plt.tight_layout()  
#                 
#     ## different type of FDR
#     
#     # No. of Observed A-neg pair / No. of expected A-neg pair
#     
#     Key_List_sort = list(np.array(Key_List)[IDX]) 
#     Pval_List_sort = list(np.array(Pval_List)[IDX]) 
#     
#     NoOfObsNeg = []
#     NoOfExpNeg = []
#     Count = 0
#     
#     for i,key in enumerate(Key_List_sort):
#         if negId1 in key or negId2 in key:
#             Count+=1
#         NoOfObsNeg.append(Count)
#         NoOfExpNeg.append(TotalNegNum/float(len(Key_List))*i)
#         
#     YY = np.array(NoOfObsNeg) / np.array(NoOfExpNeg)     
# 
#     plt.subplot(2,3,4) 
#     plt.plot(Pval_List_sort,list(YY),'r.',markersize=2,rasterized=True) 
#     plt.xlabel('P values')
#     plt.ylabel('Observed No. of neg-ctrls / Expected No. of neg-ctrls')
#     
#     plt.subplot(2,3,5)
#     IDX = np.argsort(Qval_List)
#     NoOfSig = []
#     NoOfFP = []
#     for qval in np.array(Qval_List)[IDX]:
#         NoOfSig.append(np.sum(np.array(Qval_List)<=qval))
#         NoOfFP.append(np.sum(np.array(Qval_List)<=qval) * qval)
# 
#     plt.plot(np.array(Qval_List)[IDX],NoOfSig,'r.',markersize=2,rasterized=True) 
#     plt.xlabel('Q values')
#     plt.ylabel('No. of Significant pairs')
#     plt.xlim((0,max(np.array(Qval_List)[IDX])))
#     
#     plt.subplot(2,3,6)    
#     plt.plot(NoOfSig,NoOfFP,'r.',markersize=2,rasterized=True) 
#     plt.ylabel('No. of False Positive pairs')
#     plt.xlabel('No. of Significant pairs')
#     plt.xlim((0,max(NoOfSig)))
#                         
#     plt.tight_layout()
#     plt.savefig(outbase+'_FDR_plots.pdf',format='PDF',transparent=True,dpi=1200)                    
#                                     
#     return QvalDict


def calculatePhe(enrich_dic,rec_dic,arguments):

    # Pulls in gene info, IDs, as well as GO terms (David Morgen's code)
    
    print('Retrieving gene information')
    
    # Retrieves gene information
    geneID2Name, geneID2Info, geneName2ID, geneEns2Name, geneName2Ens = retrieveInfo(mouse=arguments.mouse)
    
    # Retrieves GO information
    geneID2Comp, geneID2Proc, geneID2Fun = retrieveGO()

    if rec_dic['guide_type']=='single':
        gene2enrichment = {}
        gene2elementname = {}
        
        for element,enrichment in enrich_dic.iteritems():        
            genename = element.split(rec_dic['split_single'])[int(rec_dic['gene_idx'])]
            
            if rec_dic['neg_ctrl'] in element.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]:
                if rec_dic['neg_ctrl'] in gene2enrichment:
                    gene2enrichment[rec_dic['neg_ctrl']].append(float(enrichment))
                    gene2elementname[rec_dic['neg_ctrl']].append(element)
                else:
                    gene2enrichment[rec_dic['neg_ctrl']]=[float(enrichment)]
                    gene2elementname[rec_dic['neg_ctrl']]=[element]
            else:
                if genename in gene2enrichment:
                    gene2enrichment[genename].append(float(enrichment))
                    gene2elementname[genename].append(element)
                else:
                    gene2enrichment[genename]=[float(enrichment)]
                    gene2elementname[genename]=[element]
                    
        element_no = []
        var_all = []
        for key,value in gene2enrichment.iteritems():
            if key!=rec_dic['neg_ctrl']:
                element_no.append(len(value))
                var_all.append(np.var(value))
        
        Nctrl = np.median(element_no)
        min_var = np.percentile(var_all,5)
        
        Vctrl = np.var(gene2enrichment[rec_dic['neg_ctrl']])
        Uctrl = np.median(gene2enrichment[rec_dic['neg_ctrl']])
        ctrl_val = gene2enrichment[rec_dic['neg_ctrl']]

        # plotting histogram of control and the rest
        plt.figure(figsize=(10,10))
        n, bins, patches = plt.hist(enrich_dic.values(), 100, facecolor='k', alpha=0.2, label='all_elements')
        plt.hist(gene2enrichment[rec_dic['neg_ctrl']], bins, facecolor='m', alpha=0.5, label='neg_elements')
        plt.xlabel('element phenotypes (pZ)')
        plt.ylabel('no. of elements')
        plt.legend(loc='best')
        file_in = arguments.file_head + '_negctrl_distribution.pdf'
        plt.savefig(file_in,format='PDF',dpi=400,transparent=True)
        plt.close()
                                
        gene2no = {}
        gene2tscore = {}
        gene2medphe = {}
        gene2bootstrap = {}
        Gene2Mwu = {}
            
        count=0            
        for key,value in gene2enrichment.iteritems():
            if count==0:
                print '\nProgress',
            if count%200==0:
                print str(int(100*count/len(gene2enrichment)))+'%',
            count+=1
                            
            if key!=rec_dic['neg_ctrl']:
                Svar = (max(min_var,np.var(value)) * (len(value)-1) + Vctrl * (Nctrl-1))/(len(value)+Nctrl-2)
                ts = (np.median(value)-Uctrl)/((Svar/len(value)+Svar/Nctrl)**0.5)
                gene2no[key]=len(value)
                gene2tscore[key]=ts
                gene2medphe[key]=np.median(value)
                    
                if arguments.statistics=='all' or arguments.statistics=='bootstrap':
                    median_diff = np.median(value)-Uctrl
                    median_all = np.median(value+ctrl_val)
                    
                    value_shifted = np.array(value)-np.median(value)+median_all
                    ctrl_shifted = np.array(ctrl_val)-Uctrl+median_all
                    
                    value_bt = bootstrap_test(value_shifted,np.median,arguments.bootstrap_repeat)
                    ctrl_bt = bootstrap_test(ctrl_shifted,np.median,arguments.bootstrap_repeat)
                    
                    median_diff_bt = value_bt-ctrl_bt
                    
                    pvalue = 1-np.sum(median_diff_bt>=median_diff)/len(median_diff_bt)
                    gene2bootstrap[key]=pvalue
                    
                if arguments.statistics=='all' or arguments.statistics=='mwu':
                    pvalue = ms.mannwhitneyu(value,ctrl_val)[1]    
                    Gene2Mwu[key]=pvalue
                    
        # if arguments.statistics=='all' or arguments.statistics=='mwu':
        #     mwu_qvaldict = CalculateQvalFromPval(Gene2Mwu,rec_dic['neg_ctrl'],rec_dic['neg_ctrl'],arguments.file_head,FDR=0.05)            
        # else:
        #     mwu_qvaldict = {}                                                          

        if arguments.statistics=='all' or arguments.statistics=='mwu':
            ppval = []
            kkey = []
            for key, value in Gene2Mwu.iteritems():
                ppval.append(value)
                kkey.append(key)
                        
            adjusted_pval = multipletests(ppval,method='fdr_bh')[1]
            mwu_qvaldict={}
            for i,key in enumerate(kkey):
                mwu_qvaldict[key]=adjusted_pval[i]          
        else:
            mwu_qvaldict = {}      
                                                                                                                                                                                                                                                                                                                                                                                                                  
        ### save data (optional: save orientation specific data?)
        
        file_in = arguments.file_head + '_phenotype.csv'            
        with open(file_in, 'wb') as out_open:
            out_csv = csv.writer(out_open, delimiter=',')
            out_csv.writerow(['EnsemblID','GeneID','Symbol','GeneInfo','Localization','Process','Function','#Element','Med Phenotype (pZ)','Tscore','Bootstrap P-val ('+np.str(arguments.bootstrap_repeat)+')','MWU P-val','MWU Q-val (fdr_bh)'])
            for key, value in gene2no.iteritems():
                
                if key in geneName2Ens:
                    ensid = geneName2Ens[key]
                else:
                    ensid = 'N/A'
                if key in geneName2ID:
                    geneid = geneName2ID[key]
                else:
                    geneid = 'N/A'                
                if geneid in geneID2Info:    
                    info = geneID2Info[geneid]
                else:
                    info = 'N/A'
                if geneid in geneID2Comp:    
                    loc = geneID2Comp[geneid]
                else:
                    loc = 'N/A'
                if geneid in geneID2Proc:
                    proc = geneID2Proc[geneid]
                else:
                    proc = 'N/A'
                if geneid in geneID2Fun:
                    fun = geneID2Fun[geneid]
                else:
                    fun = 'N/A'
                if key in gene2no:
                    num = gene2no[key]
                else:
                    num = 'N/A'
                if key in gene2medphe:
                    medphe = gene2medphe[key]
                else:
                    medphe = 'N/A'
                if key in gene2tscore:
                    ts = gene2tscore[key]
                else:
                    ts = 'N/A'                        
                if key in gene2bootstrap:
                    pval = gene2bootstrap[key]
                else:
                    pval = 'N/A'                        

                if key in Gene2Mwu:
                    mwu_pval = Gene2Mwu[key]
                else:
                    mwu_pval = 'N/A' 

                if key in mwu_qvaldict:
                    mwu_qval = mwu_qvaldict[key]
                else:
                    mwu_qval = 'N/A'                                 
                                                                                                
                out_csv.writerow([ensid,geneid,key,info,loc,proc,fun,num,medphe,ts,pval,mwu_pval,mwu_qval])
        
        return gene2enrichment, gene2elementname, gene2medphe, gene2tscore, gene2bootstrap, Gene2Mwu
        
    if rec_dic['guide_type']=='double':    
  
        ### 1) calculate single knockout phenotypes of sgRNAs
        front_e_singleko = {}
        back_e_singleko = {}
    
        for key,value in enrich_dic.iteritems():
            front_e = key.split(rec_dic['split_double'])[0]
            back_e = key.split(rec_dic['split_double'])[1]
            
            if back_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]==rec_dic['neg_ctrl']:
                if front_e in front_e_singleko:
                    front_e_singleko[front_e].append(float(value))
                else:
                    front_e_singleko[front_e]=[float(value)]
                                
            if front_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]==rec_dic['neg_ctrl']:
                if back_e in back_e_singleko:
                    back_e_singleko[back_e].append(float(value))
                else:
                    back_e_singleko[back_e]=[float(value)]
    
        front_e_singleko_medphe = {}
        back_e_singleko_medphe = {}
        
        front_no = []
        for key, value in front_e_singleko.iteritems():
            front_no.append(len(value))
            
        back_no = []
        for key, value in back_e_singleko.iteritems():
            back_no.append(len(value))
        
        for key, value in front_e_singleko.iteritems():
            if len(value) >= arguments.percentile_cut * np.median(front_no): #remove element if no. of safes-element is less than percentile cut
                front_e_singleko_medphe[key]=np.median(value)
        
        for key, value in back_e_singleko.iteritems():
            if len(value) >= arguments.percentile_cut * np.median(back_no):        
                back_e_singleko_medphe[key]=np.median(value)        
                
        # calculate expected phenotypes of double-sgRNAs
        exp_phe_dic = {}
        obs_phe_dic = {}
        
        for key,value in enrich_dic.iteritems():
            front_e = key.split(rec_dic['split_double'])[0]
            back_e = key.split(rec_dic['split_double'])[1]
            
            if front_e in front_e_singleko_medphe and back_e in back_e_singleko_medphe:
                exp_phe_dic[key] = front_e_singleko_medphe[front_e] + back_e_singleko_medphe[back_e]
                obs_phe_dic[key] = value       
        
        ### 2) fitting the models        
        double_name = []; xval = []; yval = []; neg_name = []; xneg = []; yneg = [];
        
        for key in exp_phe_dic:
            double_name.append(key)        
            xval.append(exp_phe_dic[key])
            yval.append(obs_phe_dic[key])
            if rec_dic['neg_ctrl'] in key:
                xneg.append(exp_phe_dic[key])
                yneg.append(obs_phe_dic[key])
                neg_name.append(key)
    
        plt.figure(figsize=(15,15))        
        plt.subplot(2,2,1)
        plt.plot(xval,yval,'k.',label='double-element',rasterized=True,alpha=0.5)
        plt.xlabel('Expected double-element phe (pZ)')
        plt.ylabel('Observed double-element phe (pZ)')

        if arguments.fitting=='all':       
            # calculate median        
            Interval = (max(xval) - min(xval)) / 100.0
            Edge1 = np.arange(min(xval),max(xval),Interval)
            Edge0 = np.array(Edge1[0] - Interval)
            Edge2 = np.array(Edge1[-1] + Interval)
            Edge3 = np.array(Edge1[-1] + 2*Interval)
            Edge = np.hstack((Edge0,Edge1,Edge2,Edge3))
        
            Edge = Edge+Interval/2.0
            Edge = Edge[:-1]    
            XX = (Edge + Interval/2)[:-1]    
            Prct_50 = st.binned_statistic(np.array(xval), np.array(yval), statistic=lambda y: np.percentile(y,50), bins=Edge)[0]    
                
            XX = XX[~np.isnan(Prct_50)]
            YY=Prct_50[~np.isnan(Prct_50)]
            plt.plot(XX,YY,'ro',label='median')
            
            # smooth median    
            Ysmooth = running_mean(YY,5)
            Ysmooth = np.hstack((YY[0:2],Ysmooth))
            Ysmooth = np.hstack((Ysmooth,YY[-2:]))
            plt.plot(XX,Ysmooth,'y-',linewidth=1,label='smooth median')
            
            Buff = (max(xval+yval)-min(xval+yval))/20.0
        
            plt.plot([min(xval+yval)-Buff,max(xval+yval)+Buff],[min(xval+yval)-Buff,max(xval+yval)+Buff],'c-.',linewidth=3)
            plt.xlim((min(xval+yval)-Buff,max(xval+yval)+Buff))
            plt.ylim((min(xval+yval)-Buff,max(xval+yval)+Buff))
            plt.legend(loc=2)
        
            # interpolate function
            XX[-1] = max(XX[-1],max(xval))
            XX[0] = min(XX[0],min(xval))     
            f = interp1d(XX,Ysmooth)
            
        elif arguments.fitting=='neg_ctrl':       
            # calculate median
            new_xval = xval + xneg * int(len(xval)/len(xneg) * 100) # Weighing neg controls 100 times more than other data points to generate expected phenotype line
            new_yval = yval + yneg * int(len(yval)/len(yneg) * 100) # Weighing neg controls 100 times more than other data points to generate expected phenotype line
                                   
            Interval = (max(new_xval) - min(new_xval)) / 100.0
            Edge1 = np.arange(min(new_xval),max(new_xval),Interval)
            Edge0 = np.array(Edge1[0] - Interval)
            Edge2 = np.array(Edge1[-1] + Interval)
            Edge3 = np.array(Edge1[-1] + 2*Interval)
            Edge = np.hstack((Edge0,Edge1,Edge2,Edge3))
        
            Edge = Edge+Interval/2.0
            Edge = Edge[:-1]    
            XX = (Edge + Interval/2)[:-1]    
            Prct_50 = st.binned_statistic(np.array(new_xval), np.array(new_yval), statistic=lambda y: np.percentile(y,50), bins=Edge)[0]    
                
            XX = XX[~np.isnan(Prct_50)]
            YY=Prct_50[~np.isnan(Prct_50)]
            plt.plot(XX,YY,'ro',label='median')
            
            # smooth median    
            Ysmooth = running_mean(YY,5)
            Ysmooth = np.hstack((YY[0:2],Ysmooth))
            Ysmooth = np.hstack((Ysmooth,YY[-2:]))
            plt.plot(XX,Ysmooth,'y-',linewidth=1,label='smooth median')
            
            Buff = (max(xval+yval)-min(xval+yval))/20.0
        
            plt.plot([min(xval+yval)-Buff,max(xval+yval)+Buff],[min(xval+yval)-Buff,max(xval+yval)+Buff],'c-.',linewidth=3)
            plt.xlim((min(xval+yval)-Buff,max(xval+yval)+Buff))
            plt.ylim((min(xval+yval)-Buff,max(xval+yval)+Buff))
            plt.legend(loc=2)
        
            # interpolate function
            XX[-1] = max(XX[-1],max(xval))
            XX[0] = min(XX[0],min(xval))     
            f = interp1d(XX,Ysmooth)

        # calculate rawGI (sign of GI will be corrected at the gene level)
        plt.subplot(2,2,2)
        Yexp = f(xval) # Yexp is corrected expected phenotype
        temp_rawGI = np.array(yval) - Yexp # deviation
        
        if arguments.gi_option=='corr_sign':
        ## correct sign of temp_rawGI 
            rawGI = []
            for i,r_gi in enumerate(temp_rawGI):
                if xval[i]<0:
                    sign=1
                else:
                    sign=-1
                
                rawGI.append(r_gi*sign)        
        elif arguments.gi_option=='unch_sign':
            rawGI = []
            for i,r_gi in enumerate(temp_rawGI):
                sign=1
                rawGI.append(r_gi*sign)
                        
        plt.plot(xval,rawGI,'k.',rasterized=True,alpha=0.2)
        plt.plot([min(xval)-Buff,max(xval)+Buff],[0,0],'c-.',linewidth=3)
        plt.xlim((min(xval)-Buff,max(xval)+Buff))
        plt.ylim(-max(abs(np.array(rawGI)))*3,max(abs(np.array(rawGI))*3))
        plt.xlabel('expected double-element phe (pZ)')
        plt.ylabel('rawGI (pZ)')
        
        # calculate std curve of rawGI    
        Edge = [] 
        for i, xt in enumerate(sorted(xval)):
            if i%1000 == 0:
                Edge.append(xt)
        Edge.append(xt)
        
        IDX = np.digitize(xval,Edge)
        
        XGroup = {}
        YGroup = {}
    
        print '\n\nNormalizing GI........'        
        for j in range(max(IDX)):
            if j==0:
                temp_xval = np.hstack((np.array(xval)[IDX==j],np.array(xval)[IDX==j+1]))
                temp_yval = np.hstack((np.array(rawGI)[IDX==j],np.array(rawGI)[IDX==j+1]))
                XGroup[j]=temp_xval
                YGroup[j]=temp_yval            
                continue
                
            if j==max(IDX)-1:
                temp_xval = np.hstack((np.array(xval)[IDX==j],np.array(xval)[IDX==j-1]))
                temp_yval = np.hstack((np.array(rawGI)[IDX==j],np.array(rawGI)[IDX==j-1]))
                XGroup[j]=temp_xval
                YGroup[j]=temp_yval 
                continue
                
            if j>0 and j<max(IDX)-1:
                temp_xval = np.hstack((np.array(xval)[IDX==j],np.array(xval)[IDX==j-1],np.array(xval)[IDX==j+1]))
                temp_yval = np.hstack((np.array(rawGI)[IDX==j],np.array(rawGI)[IDX==j-1],np.array(rawGI)[IDX==j+1]))
                XGroup[j]=temp_xval
                YGroup[j]=temp_yval 
        
        Ystd = []
        
        for i, xt in enumerate(IDX-1):
            temp_diff = abs(XGroup[xt]-xval[i])
            sort_inds = np.argsort(temp_diff)
            Yneighbor = YGroup[xt][sort_inds][:200]
            Ystd_temp = np.std(Yneighbor)
            Ystd.append(Ystd_temp)
            
        plt.plot(XGroup[xt][sort_inds][:200],Yneighbor,'r.')
    
        # plotting std and smoothened std    
        plt.subplot(2,2,3)
        plt.plot(xval,Ystd,'k.',markersize=2,rasterized=True)

        # calculate median        
        Interval = (max(xval) - min(xval)) / 100.0
        Edge1 = np.arange(min(xval),max(xval),Interval)
        Edge0 = np.array(Edge1[0] - Interval)
        Edge2 = np.array(Edge1[-1] + Interval)
        Edge3 = np.array(Edge1[-1] + 2*Interval)
        Edge = np.hstack((Edge0,Edge1,Edge2,Edge3))
    
        Edge = Edge+Interval/2.0
        Edge = Edge[:-1]    
        XX = (Edge + Interval/2)[:-1]
        
        Prct_50 = st.binned_statistic(np.array(xval), np.array(Ystd), statistic=lambda y: np.percentile(y,50), bins=Edge)[0]    
            
        XX = XX[~np.isnan(Prct_50)]
        YY=Prct_50[~np.isnan(Prct_50)]
        plt.plot(XX,YY,'ro',label='median')    
        # smooth median    
        Ysmooth = running_mean(YY,5)
        Ysmooth = np.hstack((YY[0:2],Ysmooth))
        Ysmooth = np.hstack((Ysmooth,YY[-2:]))
        plt.plot(XX,Ysmooth,'y-',linewidth=1,label='smooth median')
        XX[-1] = max(XX[-1],max(xval))
        XX[0] = min(XX[0],min(xval))     
        f = interp1d(XX,Ysmooth)    
        Ystd_smooth = f(xval)
                        
        plt.xlim((min(xval)-Buff,max(xval)+Buff))
        plt.xlabel('expected double-element phe (pZ)')
        plt.ylabel('std calculated from rawGI (neighbor200)')
    
        # plotting normGI
        normGI = []
        for i,rgi in enumerate(rawGI):
            normGI.append(rgi/Ystd_smooth[i])
            #normGI.append(rgi/Ystd[i])
        
        plt.subplot(2,2,4)                
        plt.plot(xval,normGI,'k.',rasterized=True,alpha=0.2)
        plt.plot([min(xval)-Buff,max(xval)+Buff],[0,0],'c-.',linewidth=3)
        plt.xlim((min(xval)-Buff,max(xval)+Buff))
        plt.ylim(-max(abs(np.array(normGI)))*2,max(abs(np.array(normGI))*2))
        plt.xlabel('expected double-element phe (pZ)')
        plt.ylabel('normGI (pZ)')       
    
        # negative control plotting
        XAneg = []
        YAneg = []
        Xneg_neg = []
        Yneg_neg = []
        for i,lbl in enumerate(double_name):
            front_e = lbl.split(rec_dic['split_double'])[0]
            back_e = lbl.split(rec_dic['split_double'])[1]        
            
            if (rec_dic['neg_ctrl'] == front_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]) + (rec_dic['neg_ctrl'] == back_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]) == 1:            
                XAneg.append(xval[i])
                YAneg.append(normGI[i])
            elif (rec_dic['neg_ctrl'] == front_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]) + (rec_dic['neg_ctrl'] == back_e.split(rec_dic['split_single'])[int(rec_dic['neg_idx'])]) == 2:
                Xneg_neg.append(xval[i])
                Yneg_neg.append(normGI[i])           
    
        plt.plot(XAneg,YAneg,'m.',label='A_neg',rasterized=True,alpha=0.2)
        plt.plot(Xneg_neg,Yneg_neg,'y.',label='neg_neg',rasterized=True,alpha=0.2)
        plt.legend(loc=2)
    
        plt.plot([min(xval)-Buff,max(xval)+Buff],[0,0],'b-',linewidth=1,alpha=0.5)
        plt.plot([0,0],[-max(abs(np.array(normGI)))*2,max(abs(np.array(normGI))*2)],'b-',linewidth=1,alpha=0.5)
        plt.tight_layout()
        file_in = arguments.file_head + '_normGI_processing.pdf'
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()
        
        ### 3) statistics
        o_genepair2expphe = {} # gene pairs (orientation matters)
        o_genepair2obsphe = {} # gene pairs (orientation matters)
        o_genepair2normgi = {} # gene pairs (orientation matters)
        
        genepair2expphe = {} # gene pairs (orientation doesn't matter)
        genepair2obsphe = {} # gene pairs (orientation doesn't matter)            
        genepair2normgi = {} # gene pairs (orientation doesn't matter)    
        
        for i,double_e in enumerate(double_name):
            front_e = double_e.split(rec_dic['split_double'])[0]
            back_e = double_e.split(rec_dic['split_double'])[1]
            
            front_gene = front_e.split(rec_dic['split_single'])[int(rec_dic['gene_idx'])]
            back_gene = back_e.split(rec_dic['split_single'])[int(rec_dic['gene_idx'])]
            
            o_genepair = front_gene + '__' + back_gene
            genepair = np.sort([front_gene,back_gene])[0] + '__' + np.sort([front_gene,back_gene])[1]
            
            if o_genepair in o_genepair2normgi:
                o_genepair2expphe[o_genepair].append(xval[i])
                o_genepair2obsphe[o_genepair].append(yval[i])
                o_genepair2normgi[o_genepair].append(normGI[i])            
            else:
                o_genepair2expphe[o_genepair]=[xval[i]]
                o_genepair2obsphe[o_genepair]=[yval[i]]            
                o_genepair2normgi[o_genepair]=[normGI[i]]
            
            if genepair in genepair2normgi:
                genepair2expphe[genepair].append(xval[i])
                genepair2obsphe[genepair].append(yval[i])            
                genepair2normgi[genepair].append(normGI[i])
            else:
                genepair2expphe[genepair]=[xval[i]]
                genepair2obsphe[genepair]=[yval[i]]             
                genepair2normgi[genepair]=[normGI[i]]
        
        # median number of elements per genepair
        o_element_no = []
        o_var_all = []
        o_neg_values = []
        for key,value in o_genepair2normgi.iteritems():
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if front_g==rec_dic['neg_ctrl'] or back_g==rec_dic['neg_ctrl']:
                o_neg_values.extend(value)        
            else:        
                o_element_no.append(len(value))
                o_var_all.append(np.var(value))
        
        o_Nctrl = np.median(o_element_no)
        o_Uctrl = np.median(o_neg_values)
        o_Vctrl = np.var(o_neg_values)
        o_minVar = np.percentile(o_var_all,5)    
            
        element_no = []
        var_all = []    
        neg_values = []
        for key,value in genepair2normgi.iteritems():
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if front_g==rec_dic['neg_ctrl'] or back_g==rec_dic['neg_ctrl']:
                neg_values.extend(value)        
            else:        
                element_no.append(len(value))
                var_all.append(np.var(value))
    
        Nctrl = np.median(element_no)
        Uctrl = np.median(neg_values)
        Vctrl = np.var(neg_values)
        minVar = np.percentile(var_all,5)
                                        
        # print histograms for no. of combinations detected
        plt.figure(figsize=(15,7.5))
        plt.subplot(1,2,1)
        plt.hist(o_element_no,bins=range(0,max(o_element_no)+1),facecolor='g',alpha=0.5)
        plt.xlabel('no. of combinations per gene pair')
        plt.ylabel('no. of gene pairs')
        plt.title('orientation matters')
    
        plt.subplot(1,2,2)
        plt.hist(element_no,bins=range(0,max(element_no)+1),facecolor='g',alpha=0.5)
        plt.xlabel('no. of combinations per gene pair')
        plt.ylabel('no. of gene pairs')
        plt.title('orientation does not count')
        plt.tight_layout()
        file_in = arguments.file_head + '_histogram_no_combinations_per_genepair.pdf'
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)    
        plt.close()
        
        #get tscores and other values
    
        gene2bootstrap = {}
        o_gene2bootstrap = {}
        genepair2tscore = {}
        o_genepair2tscore = {}
        GenePair2Mwu = {}
        o_GenePair2Mwu = {}
                
        o_genepair2medexpphe = {}
        o_genepair2medobsphe = {}
        o_genepair2mednormgi = {}
                            
        for key,value in o_genepair2normgi.iteritems():
            o_Uexp = np.median(value)
            o_Nexp = min([len(value),max(o_element_no)])
            o_Vexp = max(np.var(value),o_minVar)            
            o_Svar = (o_Vexp*(o_Nexp-1)+o_Vctrl*(o_Nctrl-1))/(o_Nexp+o_Nctrl-2)            
            tscore = (o_Uexp-o_Uctrl)/((o_Svar/o_Nexp+o_Svar/o_Nctrl)**0.5)
            
            o_genepair2medexpphe[key]=np.median(o_genepair2expphe[key])
            o_genepair2medobsphe[key]=np.median(o_genepair2obsphe[key])
            o_genepair2mednormgi[key]=np.median(o_genepair2normgi[key])
            
            o_genepair2tscore[key] = tscore

        genepair2medexpphe = {}
        genepair2medobsphe = {}
        genepair2mednormgi = {}
                            
        for key,value in genepair2normgi.iteritems():
            Uexp = np.median(value)
            Nexp = min([len(value),max(element_no)])
            Vexp = max(np.var(value),minVar)            
            Svar = (Vexp*(Nexp-1)+Vctrl*(Nctrl-1))/(Nexp+Nctrl-2)            
            tscore = (Uexp-Uctrl)/((Svar/Nexp+Svar/Nctrl)**0.5)
            
            genepair2medexpphe[key]=np.median(genepair2expphe[key])
            genepair2medobsphe[key]=np.median(genepair2obsphe[key])
            genepair2mednormgi[key]=np.median(genepair2normgi[key])
            
            genepair2tscore[key] = tscore
                        
        #compare tscores of two orientations
        xx = []; yy = []; lbl1 = []; lbl2 = []; 
        for key in o_genepair2tscore:
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            flip_gp = back_g+'__'+front_g
            if flip_gp in o_genepair2tscore and (key not in lbl1+lbl2) and front_g!=back_g:
                xx.append(o_genepair2tscore[key])
                yy.append(o_genepair2tscore[flip_gp])
                lbl1.append(key)
                lbl2.append(flip_gp)
        
        Max = max(xx+yy)
        Min = min(xx+yy)    
        
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'k.',markersize=5,alpha=.5,rasterized=True)
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('tscore from AB gene pair')
        plt.ylabel('tscore from BA gene pair')
        plt.tight_layout()
        plt.legend(loc='best')            
        file_in = arguments.file_head + '_ABBA_tscore_scatterplot.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)    
        plt.close()
                            
        if arguments.statistics=='all' or arguments.statistics=='bootstrap':
            count=0
            start_time = timeit.default_timer()            
            for key,value in o_genepair2normgi.iteritems():
                if count==0:                    
                    print '\nProgress for o_bootstrap pval',
                if count%200==0:
                    elapsed = timeit.default_timer() - start_time                    
                    print str(int(100*count/len(o_genepair2normgi)))+'%', int(elapsed), 'sec',
                count+=1        
                                            
                median_diff = np.median(value)-o_Uctrl
                median_all = np.median(value+o_neg_values)
                
                value_shifted = np.array(value)-np.median(value)+median_all
                ctrl_shifted = np.array(o_neg_values)-o_Uctrl+median_all
                
                value_bt = bootstrap_test(np.float32(value_shifted),np.median,arguments.bootstrap_repeat)                
                ctrl_bt = bootstrap_test(np.float32(ctrl_shifted),np.median,arguments.bootstrap_repeat)
                
                median_diff_bt = value_bt-ctrl_bt
                
                pvalue = 1-np.sum(median_diff_bt>=median_diff)/len(median_diff_bt)
                o_gene2bootstrap[key]=pvalue         
                                        
            count=0
            start_time = timeit.default_timer()                        
            for key,value in genepair2normgi.iteritems():
                if count==0:
                    print '\nProgress for bootstrap pval',
                if count%200==0:
                    elapsed = timeit.default_timer() - start_time                                        
                    print str(int(100*count/len(genepair2normgi)))+'%', int(elapsed), 'sec',
                count+=1 
                                    
                median_diff = np.median(value)-Uctrl
                median_all = np.median(value+neg_values)
                
                value_shifted = np.array(value)-np.median(value)+median_all
                ctrl_shifted = np.array(neg_values)-Uctrl+median_all
                
                value_bt = bootstrap_test(np.float32(value_shifted),np.median,arguments.bootstrap_repeat)
                ctrl_bt = bootstrap_test(np.float32(ctrl_shifted),np.median,arguments.bootstrap_repeat)
                
                median_diff_bt = value_bt-ctrl_bt
                
                pvalue = 1-np.sum(median_diff_bt>=median_diff)/len(median_diff_bt)
                gene2bootstrap[key]=pvalue                                               
    
            #compare pvalues of two orientations
                                                                
            xx = []; yy = []; lbl1 = []; lbl2 = []; 
            for key in o_gene2bootstrap:
                front_g = key.split('__')[0]
                back_g = key.split('__')[1]
                
                flip_gp = back_g+'__'+front_g
                if flip_gp in o_gene2bootstrap and (key not in lbl1+lbl2) and front_g!=back_g:
                    value1 = o_gene2bootstrap[key]
                    value2 = o_gene2bootstrap[flip_gp]
                    
                    if value1==0:
                        temp=list(set(o_gene2bootstrap.values()))
                        value1=np.sort(temp)[1]*0.1
                    if value2==0:
                        temp=list(set(o_gene2bootstrap.values()))
                        value2=np.sort(temp)[1]*0.1
                    
                    xx.append(np.log10(value1)*-1*np.sign(o_genepair2tscore[key]))
                    yy.append(np.log10(value2)*-1*np.sign(o_genepair2tscore[flip_gp]))
                    lbl1.append(key)
                    lbl2.append(flip_gp)
            
            Max = max(xx+yy)
            Min = min(xx+yy)    
            
            Slope = np.polyfit(xx,yy,1)[0]
            Intercept = np.polyfit(xx,yy,1)[1]
            
            plt.figure(figsize=(10,10))
            plt.plot(xx,yy,'k.',markersize=5,alpha=.2,rasterized=True)
            Buff= (Max-Min)/5.0
            plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
            plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
            plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
            plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
            
            PearsonCoeff=st.stats.pearsonr(xx,yy)
            plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
            
            plt.xlim([Min-Buff,Max+Buff])
            plt.ylim([Min-Buff,Max+Buff])
            plt.xlabel('signed log10(pval) from AB gene pair')
            plt.ylabel('signed log10(pval) from BA gene pair')
            plt.legend(loc='best')
            plt.tight_layout()
            file_in = arguments.file_head + '_ABBA_bootstrap_signed_log10pval_scatterplot.pdf'        
            plt.savefig(file_in,format='PDF',transparent=True,dpi=400)                                                                        
            plt.close()

        if arguments.statistics=='all' or arguments.statistics=='mwu':
            count=0
            start_time = timeit.default_timer()            
            for key,value in o_genepair2normgi.iteritems():
                if count==0:                    
                    print '\nProgress for o_mwu pval',
                if count%200==0:
                    elapsed = timeit.default_timer() - start_time                    
                    print str(int(100*count/len(o_genepair2normgi)))+'%', int(elapsed), 'sec',
                count+=1        
                
                pvalue = ms.mannwhitneyu(value,o_neg_values)[1]    
                o_GenePair2Mwu[key]=pvalue         
                                        
            count=0
            start_time = timeit.default_timer()                        
            for key,value in genepair2normgi.iteritems():
                if count==0:
                    print '\nProgress for o_mwu pval',
                if count%200==0:
                    elapsed = timeit.default_timer() - start_time                                        
                    print str(int(100*count/len(genepair2normgi)))+'%', int(elapsed), 'sec',
                count+=1 
                                    
                pvalue = ms.mannwhitneyu(value,neg_values)[1]    
                GenePair2Mwu[key]=pvalue                                               
    
            #compare pvalues of two orientations                                                                
            xx = []; yy = []; lbl1 = []; lbl2 = []; 
            for key in o_GenePair2Mwu:
                front_g = key.split('__')[0]
                back_g = key.split('__')[1]
                
                flip_gp = back_g+'__'+front_g
                if flip_gp in o_GenePair2Mwu and (key not in lbl1+lbl2) and front_g!=back_g:
                    value1 = o_GenePair2Mwu[key]
                    value2 = o_GenePair2Mwu[flip_gp]
                    
                    if value1==0:
                        temp=list(set(o_GenePair2Mwu.values()))
                        value1=np.sort(temp)[1]*0.1
                    if value2==0:
                        temp=list(set(o_GenePair2Mwu.values()))
                        value2=np.sort(temp)[1]*0.1
                    
                    xx.append(np.log10(value1)*-1*np.sign(o_genepair2tscore[key]))
                    yy.append(np.log10(value2)*-1*np.sign(o_genepair2tscore[flip_gp]))
                    lbl1.append(key)
                    lbl2.append(flip_gp)
            
            Max = max(xx+yy)
            Min = min(xx+yy)    
            
            Slope = np.polyfit(xx,yy,1)[0]
            Intercept = np.polyfit(xx,yy,1)[1]
            
            plt.figure(figsize=(10,10))
            plt.plot(xx,yy,'k.',markersize=5,alpha=.2,rasterized=True)
            Buff= (Max-Min)/5.0
            plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
            plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
            plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
            plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
            
            PearsonCoeff=st.stats.pearsonr(xx,yy)
            plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
            
            plt.xlim([Min-Buff,Max+Buff])
            plt.ylim([Min-Buff,Max+Buff])
            plt.xlabel('signed log10(pval) from AB gene pair')
            plt.ylabel('signed log10(pval) from BA gene pair')
            plt.legend(loc='best')
            plt.tight_layout()
            file_in = arguments.file_head + '_ABBA_mwu_signed_log10pval_scatterplot.pdf'        
            plt.savefig(file_in,format='PDF',transparent=True,dpi=400)                                                                        
            plt.close()
        
        ### 4) save data (optional: save orientation specific data?)
        # if arguments.statistics=='all' or arguments.statistics=='mwu':
        #     mwu_qvaldict = CalculateQvalFromPval(GenePair2Mwu,rec_dic['neg_ctrl'],rec_dic['neg_ctrl'],arguments.file_head,FDR=0.05)
        # else:
        #     mwu_qvaldict = {}        
        if arguments.statistics=='all' or arguments.statistics=='mwu':
            ppval = []
            kkey = []
            for key, value in GenePair2Mwu.iteritems():
                ppval.append(value)
                kkey.append(key)
                        
            adjusted_pval = multipletests(ppval,method='fdr_bh')[1]
            mwu_qvaldict={}
            for i,key in enumerate(kkey):
                mwu_qvaldict[key]=adjusted_pval[i]          
        else:
            mwu_qvaldict = {}
        
        # 01.14.2020 update !!!: p values and q values are added for double knockout phenotypes for Attardi lab. 
        double_neg = rec_dic['neg_ctrl'] + '__' + rec_dic['neg_ctrl']
        neg_phe_values = genepair2obsphe[double_neg]        
        
        count = 0
        GenePair2Mwu_Phe = {}        
        for key,value in genepair2obsphe.iteritems():
            if count==0:
                print '\nProgress for mwu pval on double ko phe',
            if count%200==0:
                print str(int(100*count/len(genepair2obsphe)))+'%', int(elapsed), 'sec',
            count+=1 
                                
            pvalue = ms.mannwhitneyu(value,neg_phe_values)[1]    
            GenePair2Mwu_Phe[key]=pvalue
            
        if arguments.statistics=='all' or arguments.statistics=='mwu':
            ppval = []
            kkey = []
            for key, value in GenePair2Mwu_Phe.iteritems():
                ppval.append(value)
                kkey.append(key)
                        
            adjusted_pval = multipletests(ppval,method='fdr_bh')[1]
            mwu_qvaldict_phe={}
            for i,key in enumerate(kkey):
                mwu_qvaldict_phe[key]=adjusted_pval[i]          
        else:
            mwu_qvaldict_phe = {}
                                           
        file_in = arguments.file_head + '_phenotype.csv'
        with open(file_in, 'wb') as out_open:
            out_csv = csv.writer(out_open, delimiter=',')
            out_csv.writerow(['Gene_pair','GeneA Info','GeneB Info','GeneA Process','GeneB Process','GeneA single phe (pZ)','No. elements (GeneA)','GeneB single phe (pZ)','No. elements (GeneB)',\
            'GeneAB double exp phe (pZ)', 'GeneAB double obs phe (pZ)','No. of elements for GeneAB double','MWU Pval for phe','MWU Qval for phe (fdr_bh)',\
            'NormGI (sign not corrected)', 'Tscore (sign corrected)', 'Bootstrap Pval', 'MWU Pval','MWU Qval (fdr_bh)'])
            
            for key in genepair2tscore:
                front_g = key.split('__')[0]
                back_g = key.split('__')[1]
    
                if front_g in geneName2ID:
                    f_geneid = geneName2ID[front_g]
                else:
                    f_geneid = 'N/A'
                if back_g in geneName2ID:
                    b_geneid = geneName2ID[back_g]
                else:
                    b_geneid = 'N/A'                
                            
                if f_geneid in geneID2Info:    
                    f_info = geneID2Info[f_geneid]
                else:
                    f_info = 'N/A'
                if b_geneid in geneID2Info:    
                    b_info = geneID2Info[b_geneid]
                else:
                    b_info = 'N/A'
                    
                if f_geneid in geneID2Proc:    
                    f_proc = geneID2Proc[f_geneid]
                else:
                    f_proc = 'N/A'
                if b_geneid in geneID2Proc:    
                    b_proc = geneID2Proc[b_geneid]
                else:
                    b_proc = 'N/A'                
                    
                front_single = np.sort([front_g,rec_dic['neg_ctrl']])[0]+'__'+np.sort([front_g,rec_dic['neg_ctrl']])[1]
                if front_single in genepair2medobsphe:
                    front_sko = genepair2medobsphe[front_single]
                    front_no = len(genepair2obsphe[front_single])
                else:
                    front_sko = 'N/A'
                    front_no = 'N/A'
                    
                back_single = np.sort([back_g,rec_dic['neg_ctrl']])[0]+'__'+np.sort([back_g,rec_dic['neg_ctrl']])[1]
                if back_single in genepair2medobsphe:
                    back_sko = genepair2medobsphe[back_single]
                    back_no = len(genepair2obsphe[back_single])
                else:
                    back_sko = 'N/A'
                    back_no = 'N/A'
                    
                doubleexp = genepair2medexpphe[key]
                doubleobs = genepair2medobsphe[key]
                double_no = len(genepair2expphe[key])
                doublephe_mwu = GenePair2Mwu_Phe[key]
                doublephe_qval = mwu_qvaldict_phe[key]
    
                norm_gi = genepair2mednormgi[key]
                if key in genepair2tscore:
                    t_score = genepair2tscore[key]
                else:
                    t_score = 'N/A'
                if key in gene2bootstrap:    
                    boot_strap = gene2bootstrap[key]
                else:
                    boot_strap = 'N/A'
                if key in GenePair2Mwu:    
                    mwu = GenePair2Mwu[key]
                else:
                    mwu = 'N/A'
                if key in mwu_qvaldict:    
                    mwu_q = mwu_qvaldict[key]
                else:
                    mwu_q = 'N/A'        
                
                out_csv.writerow([key,f_info,b_info,f_proc,b_proc,front_sko,front_no,back_sko,back_no,doubleexp,doubleobs,double_no,doublephe_mwu,doublephe_qval,\
                norm_gi,t_score,boot_strap,mwu,mwu_q])
                                                                                                
        ### 5) plotting single knockout phenotypes of A_safe, safe_A genepairs
        xx=[]; yy=[];
                
        for key,value in o_genepair2medobsphe.iteritems():
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if back_g == rec_dic['neg_ctrl'] and back_g+'__'+front_g in o_genepair2medobsphe:
                xx.append(value)
                yy.append(o_genepair2medobsphe[back_g+'__'+front_g])
        
        Max = max(xx+yy)
        Min = min(xx+yy)    
        
        Slope = np.polyfit(xx,yy,1)[0]
        Intercept = np.polyfit(xx,yy,1)[1]
        
        plt.figure(figsize=(10,10))
        plt.plot(xx,yy,'k.',markersize=5,alpha=1,rasterized=True)
        Buff= (Max-Min)/5.0
        plt.plot([Min-Buff,Max+Buff],[0,0],'k-',alpha=0.5)        
        plt.plot([0,0],[Min-Buff,Max+Buff],'k-',alpha=0.5)
        plt.plot([Min-Buff,Max+Buff],[Min-Buff,Max+Buff],'k-',alpha=0.25)
        plt.plot([Min-Buff,Max+Buff],[(Min-Buff)*Slope+Intercept,(Max+Buff)*Slope+Intercept],'b-.',linewidth=2,label='Y='+str(Slope)[0:4]+'X+'+str(Intercept)[0:4])
        
        PearsonCoeff=st.stats.pearsonr(xx,yy)
        plt.text(Min+(Max-Min)/5,Max-(Max-Min)/10,'PearsonCoeff='+str(round(PearsonCoeff[0],3)))
        
        plt.xlim([Min-Buff,Max+Buff])
        plt.ylim([Min-Buff,Max+Buff])
        plt.xlabel('A_safe phe (pZ)')
        plt.ylabel('safe_A phe (pZ)')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = arguments.file_head + '_A_safe_vs_safe_A_phe_scatterplot.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)                
        plt.close()
        
    ### 6) plotting double knockout phenotypes of A_B, B_A genepairs
        xx=[]; yy=[]; lbl1=[]; lbl2=[];
                
        for key,value in o_genepair2medobsphe.iteritems():
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if back_g+'__'+front_g in o_genepair2medobsphe and back_g!=front_g and key not in lbl1 and key not in lbl2:
                xx.append(value)
                yy.append(o_genepair2medobsphe[back_g+'__'+front_g])

            lbl1.append(key)
            lbl2.append(back_g+'__'+front_g) 
                
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
        plt.xlabel('A_B gene phe (pZ)')
        plt.ylabel('B_A gene phe (pZ)')
        plt.tight_layout()
        plt.legend(loc='best')
        
        file_in = arguments.file_head + '_AB_vs_BA_gene_phe_scatterplot.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()
        
    ### 7) plotting norm_gi
        xx=[]; yy=[]; lbl1=[]; lbl2=[];
                
        for key,value in o_genepair2mednormgi.iteritems():
            front_g = key.split('__')[0]
            back_g = key.split('__')[1]
            
            if back_g+'__'+front_g in o_genepair2mednormgi and back_g!=front_g and key not in lbl1 and key not in lbl2:
                xx.append(value)
                yy.append(o_genepair2mednormgi[back_g+'__'+front_g])

            lbl1.append(key)
            lbl2.append(back_g+'__'+front_g)                
        
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
        plt.xlabel('A_B norm_gi')
        plt.ylabel('B_A norm_gi')
        plt.tight_layout()
        plt.legend(loc='best')        
        file_in = arguments.file_head + '_AB_vs_BA_normgi_scatterplot.pdf'        
        plt.savefig(file_in,format='PDF',transparent=True,dpi=400)
        plt.close()
        
        return front_e_singleko, back_e_singleko, front_e_singleko_medphe, back_e_singleko_medphe,\
        double_name, xval, yval, rawGI, normGI, neg_name, xneg, yneg,\
        o_genepair2expphe, o_genepair2obsphe, o_genepair2normgi,\
        o_genepair2medexpphe, o_genepair2medobsphe, o_genepair2mednormgi,\
        genepair2expphe, genepair2obsphe, genepair2normgi,\
        genepair2medexpphe, genepair2medobsphe, genepair2mednormgi,\
        o_gene2bootstrap, o_genepair2tscore, o_GenePair2Mwu, gene2bootstrap, genepair2tscore, GenePair2Mwu
        
