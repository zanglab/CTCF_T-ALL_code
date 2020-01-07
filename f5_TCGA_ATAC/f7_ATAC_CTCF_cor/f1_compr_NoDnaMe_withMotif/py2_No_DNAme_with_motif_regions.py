
import os,re,argparse,sys
import bisect
import pandas as pd
from scipy import stats
import association_with_regions
import numpy as np
import CTCF_TALL_modules_new


      


    
def main():    
    
        
    methylated_df = pd.read_csv('Jurkat_DNA_methy_status_GT15.csv',index_col=0)    
#     methylated_df = pd.read_csv('Jurkat_DNA_methy_status_GT15_CpGs_GT3.csv',index_col=0)
    
    union_df = CTCF_TALL_modules_new.return_occupancy_filtered()
    union_df = union_df[union_df['motif_strand']!='N']
    
    union_df = union_df.loc[union_df.index.difference(methylated_df.index)]
    union_df.to_csv('union_CTCF_No_Jurkat_DNAmethylation_with_motif.csv')
    union_df.insert(3,'id',union_df.index)
    union_df.to_csv('union_CTCF_No_Jurkat_DNAmethylation_with_motif.bed',sep='\t',header=None,index=False)
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()

