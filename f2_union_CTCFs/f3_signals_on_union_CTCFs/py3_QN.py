import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from operator import itemgetter
#def expand_region(summitlist):
#from reads_count import read_count_on_mapfile

def quantile_normalization(dataframe):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = dataframe
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized

def main(indir):               
    
    df = pd.read_csv('f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings.csv',sep='\t',index_col=0)
    df = quantile_normalization(df)
    df = df.round(2)
    df.to_csv('f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized.csv',sep='\t') 
    
    
    
          
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    #parser.add_argument('-b', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input file to be compared as basic', metavar = '<file>')
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main(args.indir)