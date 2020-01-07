import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from GenomeData import *

def main():
    
    ctcf_domain_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo.csv'
    with open(ctcf_domain_file) as ctcf_domain_inf:
        df = pd.read_csv(ctcf_domain_inf,sep='\t')
    # left position of domain
    df['left100k'] = df['middle']-100000  
    df['domain_100k_left'] = df[['left100k','domain_left']].min(axis=1)
    df['domain_100k_left'] = df['domain_100k_left'].clip(0) # set negative values to 0
    df['left_1M'] = df['middle']-1000000
    df['domain_100k_1M_left'] = df[['domain_100k_left','left_1M']].max(axis=1)
    
    # right position of domain
    df['right100k'] = df['middle']+100000 
    df['domain_100k_right'] = df[['right100k','domain_right']].max(axis=1)  
    df['domain_100k_right'] = df['domain_100k_right'].clip(0)
    df['right_1M'] = df['middle']+1000000
    df['domain_100k_1M_right'] = df[['domain_100k_right','right_1M']].min(axis=1)
    
    df.to_csv('all_CTCF_domainInfo_full.csv',index=False)
    
    # only keep those domain boundaries (left|right) that are NOT modified by 100k/1M limitation
    left_positions = df[df['domain_left']==df['domain_100k_1M_left']][['chr','domain_left','id']].rename(columns = {'domain_left':'boundary'})
    right_positions = df[df['domain_right']==df['domain_100k_1M_right']][['chr','domain_right','id']].rename(columns = {'domain_right':'boundary'})
    df_boundary= pd.concat([left_positions,right_positions])
    df_boundary['boundary']=df_boundary['boundary'].astype(int)
    df_boundary.to_csv('all_CTCF_domainInfo_raw_boundaries.csv',sep=',',index=False) 

    # keep unique position and write out bed file
    df_boundary = df_boundary.drop_duplicates(subset=['chr','boundary'])
    df_boundary = df_boundary[['chr','boundary']]
    df_boundary['boundaryPlus'] = df_boundary['boundary']+1
    df_boundary.to_csv('all_CTCF_domainInfo_raw_boundaries.bed',sep='\t',index=False,header=None) 

            

 


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--geneid', action = 'store', type = str,dest = 'geneid', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
