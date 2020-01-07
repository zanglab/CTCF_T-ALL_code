import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from GenomeData import *

def main():
    
    hic_df = pd.DataFrame()
    
    hic_boundary_file = 'NIHMS828671-supplement-3.xlsx'
    xls = pd.ExcelFile(hic_boundary_file)
    for sheet_name in xls.sheet_names:
        hic_boundary_df = pd.read_excel(hic_boundary_file,sheet_name = sheet_name,header=None)
        hic_df = pd.concat([hic_df,hic_boundary_df])
    
    hic_df = hic_df.drop_duplicates()
    hic_df.to_csv('hic_boundary.bed',sep='\t',index=False,header=None)     
        
        
            

 


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
