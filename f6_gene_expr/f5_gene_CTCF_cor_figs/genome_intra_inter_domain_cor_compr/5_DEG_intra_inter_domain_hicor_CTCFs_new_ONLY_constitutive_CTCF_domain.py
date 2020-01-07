'''
This is used to show that hicor gene-CTCF tends to occur within same domain
there still high-cor gene-CTCF pairs genome wide

for each gene, 

plot the distribution of distances from gene-CTCF pair that
with R-squared values higher than a threshold
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})
import random
import CTCF_TALL_modules_new

chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
# pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/fz_matrices'

def intra_inter_cor_collection(regions,outdir):
    
    os.makedirs(outdir,exist_ok=True)
    for chrom in chroms:
        df = pd.DataFrame()
        dis_df,domain_df,cor_df = CTCF_TALL_modules_new.return_dis_domain_cor_df(chrom)
        for start_index in np.arange(len(regions)-1):
            intra_values,inter_values = np.array([]),np.array([])
            a,b = regions[start_index],regions[start_index+1]
            #print(outdir,a,b)
            intra_df = cor_df[domain_df[(dis_df>a)&(dis_df<=b)]==2]
            inter_df = cor_df[domain_df[(dis_df>a)&(dis_df<=b)]==-1]
            for column in intra_df.columns:
                intra_values = np.append(intra_values,intra_df[column].dropna().values.tolist())
                inter_values = np.append(inter_values,inter_df[column].dropna().values.tolist())
            df.loc[b,'intra_all'] = len(intra_values)
            df.loc[b,'intra_GT0.25'] = len(intra_values[np.power(intra_values,2)>0.25])
            df.loc[b,'inter_all'] = len(inter_values)
            df.loc[b,'inter_GT0.25'] = len(inter_values[np.power(inter_values,2)>0.25])
        df.to_csv(outdir+os.sep+'{}.csv'.format(chrom))
            
def main():

    
    outdir = 'intra_inter_domain_CTCF_gene_cor_compr'
    os.makedirs(outdir,exist_ok=True)
    
    regions = np.append(np.arange(0,100000,2500),np.logspace(5,7,num=201))
    #print(regions)
    outdir_new = outdir+os.sep+'sep_2p5k_201'  
    intra_inter_cor_collection(regions,outdir_new)


    regions = np.append(np.arange(0,100000,2000),np.logspace(5,7,num=301))
    #print(regions)
    outdir_new = outdir+os.sep+'sep_2k_301'  
    intra_inter_cor_collection(regions,outdir_new)



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
  
    main()
