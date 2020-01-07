import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new





def main():

    outdir = 'f2_T_ALL_gained_lost_SMC3_binding'
    os.makedirs(outdir,exist_ok=True)


    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding("T-ALL")

    infile = 'f1_binding_pattern_csv/SMC3_CUTLL1_on_union.csv'    
    with open(infile) as inf:
        binding_df = pd.read_csv(inf,index_col=0,sep='\t')
        

    tmp_df = binding_df.loc[gained_df.index]
    tmp_df.to_csv(outdir+os.sep+'CUTLL1_SMC3_on_T-ALL_gained.csv')

   
    tmp_df = binding_df.loc[lost_df.index]
    tmp_df.to_csv(outdir+os.sep+'CUTLL1_SMC3_on_T-ALL_lost.csv')
   


    
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

