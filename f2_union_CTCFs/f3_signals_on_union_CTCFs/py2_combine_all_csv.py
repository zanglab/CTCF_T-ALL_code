import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from operator import itemgetter
#def expand_region(summitlist):
#from reads_count import read_count_on_mapfile



def main(indir):               
    
    outdir='f2_signals_on_union_bindings'
    os.makedirs(outdir,exist_ok=True)
    
    indir = 'f1_RPKM_csv'
    infiles = glob.glob(indir+os.sep+'*.csv') 
    df_all = pd.DataFrame()
    for infile in infiles:
        name = os.path.basename(infile).split('.csv')[0];print(name)
        df = pd.read_csv(infile,sep='\t',header=None,index_col=0,names = ['{}'.format(name)])
        df_all = pd.concat([df_all,df],axis=1)
#     df_all = df_all.fillna(0)
    df_all = df_all.round(2)
    #print(df_all)
    df_all.to_csv('{}/signals_RPKM_on_all_CTCF_bindings.csv'.format(outdir),sep='\t') 
    
    
    
    
          
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