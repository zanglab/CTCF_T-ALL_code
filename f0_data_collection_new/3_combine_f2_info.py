import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import multiprocessing
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
    ctcf_gsm = pd.read_csv('NCBI_ALL_CTCF_GSE_GSM.csv',sep='\t')
    all_gsms = ctcf_gsm['GSM'].tolist()
    old = pd.read_excel('f0_infile/CTCF_377_cellLine_tissue_20180320_yy.xlsx',sep='\t')
    old_gsms = old['GSM'].tolist()
    test = set(old_gsms).difference(all_gsms)
    # test = {nan, 'GSM1224675', 'GSM1817661', 'GSM1224674', 'GSM970216', 'GSM1224673', 'SRR346402', 'SRR346403', 'GSM1817657', 'GSM1224672', 'SRR346401', 'GSM1138985'}
    # add the 4 datasets in GSE50611 back
    print(test);exit(0)
    '''
    
    ctcf_gsm = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_201906.csv',sep='\t',index_col=0,header=None)
    ctcf_gsm.columns=['GSE']
    # get the tissue info for each label
    exist_data = pd.read_csv('NCBI_Homo_ChIPseq_GSE_GSM_allInfo.csv.dropdup.saved',sep='\t',index_col=0)
    ctcf_gsm = ctcf_gsm.loc[ctcf_gsm.index.difference(exist_data.index)]

    i=0
    for gsm in ctcf_gsm.index:
        
        if os.path.isfile('f2_gsmInfo/{}.txt'.format(gsm)):
            with open('f2_gsmInfo/{}.txt'.format(gsm),'r') as inf:
                if i%50==0:
                    print(i)
                i+=1
                lines = inf.readlines()[0]
                #print(lines.strip('\n').split('\t'));exit()
                gsm,antibody,title,sourcename,organism,characteristics,library = lines.strip('\n').split('\t')
                ctcf_gsm.loc[gsm,'antibody'] = antibody
                ctcf_gsm.loc[gsm,'library'] = library
                ctcf_gsm.loc[gsm,'organism'] = organism
                ctcf_gsm.loc[gsm,'title'] = title
                ctcf_gsm.loc[gsm,'sourcename'] = sourcename 
                ctcf_gsm.loc[gsm,'characteristics'] = characteristics
            #print(ctcf_gsm.dropna());exit(0)
    print(ctcf_gsm)
    ctcf_gsm = pd.concat([ctcf_gsm,exist_data])
    ctcf_gsm.to_csv('NCBI_Homo_ChIPseq_GSE_GSM_allInfo.csv',sep='\t',index=True)
    ctcf_gsm = ctcf_gsm.loc[~ctcf_gsm.index.duplicated()] # drop duplicates
    ctcf_gsm.to_csv('NCBI_Homo_ChIPseq_GSE_GSM_allInfo.csv.dropdup.saved',sep='\t',index=True)
    
    

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
