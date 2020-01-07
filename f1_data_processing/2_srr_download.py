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

#     outdir='f2_GSM2SRR'
#     os.makedirs(outdir,exist_ok=True)
     
    df = pd.read_csv('f1_GSM2SRR/ctcf_new_GSM_to_SRR.csv',sep=',',index_col=0)    
    for gsmid in df.index:
#         print('##',gsmid)
        srrs = df.loc[gsmid,'SRR'].split('&&')
        layout = df.loc[gsmid,'layout']
        for SRRid in srrs:
#             print('prefetch -v {} -X 30G'.format(SRRid))   
#             print('wait\n')

## check file exist
            srr_file='/nv/vol190/zanglab/zw5j/data/ncbi/public/sra/{}.sra'.format(SRRid)
            if not os.path.exists(srr_file):
               print(gsmid, SRRid)

    
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
