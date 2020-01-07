import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
#from get_reads_positions import reads_positions
import CTCF_TALL_modules_new

def binding_signal_in_union():
    signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized.csv'
#     signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized_head10000.csv'
    with open(signal_file) as signal_inf:
        signal_df = pd.read_csv(signal_inf,sep='\t',index_col=0)
    return signal_df

   
def main():

    outdir='f1_signal_compr_layout_TALL_cases'
    os.makedirs(outdir,exist_ok=True)  
    signal_df = binding_signal_in_union()

    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL');#print(gained_df)
    tall_bindings = gained_df.index.union(lost_df.index);print(tall_bindings)
    signal_df = signal_df.loc[tall_bindings].dropna()
    signal_df.to_csv(outdir+os.sep+'signals_RPKM_QuantileNormalized.TALL_binding_filtered.csv')



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
