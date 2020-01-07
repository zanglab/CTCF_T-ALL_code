import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
# import multiprocessing
#from GenomeData import *
#import scipy
#from scipy import stats
#from scipy.cluster.hierarchy import linkage,dendrogram,cut_tree
#sys.setrecursionlimit(10000)
# import urllib.request
# from bs4 import BeautifulSoup
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

    ctcf_df = pd.read_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue_withhJournal.csv',sep='\t',index_col=0)
    
    # CTCF peaks
    peaks_df = pd.read_csv('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f3_union_CTCF_regions/f2_union_all_peaks/f1_peak2000_datasets_and_union_summits/num_peaks.csv',sep='\t',index_col=0)
    peaks_df_new = pd.read_csv('/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_CTCF_peaks.csv',index_col=0)
    peaks_df_new.columns = ['peak_nums']
    peaks_df = pd.concat([peaks_df,peaks_df_new]).dropna()#;print(peaks_df.shape)
    
#     blacklist = pd.read_csv('blacklist.log',sep='\t',index_col=0,header=None)
    blacklist = pd.read_excel('BalckList_new.xlsx',index_col=0)
    blacklist = blacklist[blacklist['status']=='deleted']#;print(blacklist)
    
    # add peaks info
    ctcf_df = ctcf_df.merge(peaks_df,left_index=True,right_index=True,how='left')
    ctcf_df = ctcf_df.loc[ctcf_df.index.difference(blacklist.index)]
    ctcf_df[ctcf_df['peak_nums']>0].to_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue_withhJournal_peaks_final.csv')
    ctcf_df[ctcf_df['peak_nums']>2000].to_csv('CTCF_ChIPseq_GSE_GSM_description_addTissue_withhJournal_peaks_final_peak2000.csv')






    
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
