import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import multiprocessing
from collections import Counter
#from GenomeData import *
#import scipy
#from scipy import stats
#from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
#sys.setrecursionlimit(10000)
import urllib.request
from bs4 import BeautifulSoup
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
                                

def main(infile):
    # this is to check if the 377 CTCF datasets in the newly collected datasets
    '''
    ctcf_gsm = pd.read_excel('CTCF_NCBI_ALLCTCF_manu_revised_tissue_IF5.xlsx',sep='\t')
    all_gsms = ctcf_gsm['GSM'].tolist();print(len(all_gsms))
    old = pd.read_excel('f0_infile/CTCF_377_cellLine_tissue_20180320_yy.xlsx',sep='\t')
    old_gsms = old['GSM'].tolist()
    test = set(old_gsms).difference(all_gsms)
    # test = {nan, 'GSM1224675', 'GSM1817661', 'GSM1224674', 'GSM970216', 'GSM1224673', 'SRR346402', 'SRR346403', 'GSM1817657', 'GSM1224672', 'SRR346401', 'GSM1138985'}
    # add the 4 datasets in GSE50611 back
    print(test);exit(0)
    
    '''
    
    
    ctcf_gsm = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_allInfo.csv',sep='\t',index_col=0)
    ctcf_gsm = ctcf_gsm.loc[~ctcf_gsm.index.duplicated()] # drop duplicates, keep one
    ctcf_gsm = ctcf_gsm[(ctcf_gsm['organism']=='Homo sapiens')]
    print('ALL Homo:\t',ctcf_gsm.shape)
    print('\n==========')
    counter = Counter(ctcf_gsm['library'])
    print(counter.most_common(10))
    print('==========\n')
    ctcf_gsm = ctcf_gsm[(ctcf_gsm['library']=='ChIP-Seq')]
    print('Homo & ChIP-seq filtered:\t',ctcf_gsm.shape)
    
    keep_gsm, drop_gsm = [],[]
    for gsm in ctcf_gsm.index:
        # keep those with CTCF in titles
        #print(gsm,ctcf_gsm.loc[gsm,'title'])
        if re.search('CTCF',ctcf_gsm.loc[gsm,'title']):
            keep_gsm.append(gsm)
        # keep those with CTCF in antibody   
        if not pd.isnull(ctcf_gsm.loc[gsm,'characteristics']):
                if re.search('CTCF',ctcf_gsm.loc[gsm,'characteristics'],flags=re.I):
                    keep_gsm.append(gsm)
        
        # delete those with ['none','NA','N/A'] in antibody
        #for empty in ['none','NA','N/A']:
            #if not pd.isnull(ctcf_gsm.loc[gsm,'antibody']):
                #if re.search(empty,ctcf_gsm.loc[gsm,'antibody'],flags=re.I):
                    #drop_gsm.append(gsm)
    #print(set(keep_gsm).intersection(drop_gsm))
    ctcf_gsm = ctcf_gsm.loc[set(keep_gsm)]
    #ctcf_gsm = ctcf_gsm.loc[keep_gsm]
    print('then CTCF filtered:\t',ctcf_gsm.shape)
    
    ctcf_gsm.to_csv('NCBI_Homo_ChIPseq_GSE_GSM_CTCF_filtered.csv',sep='\t')
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)
