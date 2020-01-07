'''
This file is used to compr the CTCF binding signal 
in cancer cell lines
between gained/lost/ctrl for each cancer type
'''

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

    outdir='f3_gained_constitutive_CTCF_binding_signal_compr'
    os.makedirs(outdir,exist_ok=True)  
    
    signal_df = binding_signal_in_union()
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()
#     print(signal_df);exit()
     
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype);#print(gained_df)
        cancercols,normalcols = CTCF_TALL_modules_new.cancer_specific_cancer_normal_gsmID(cancertype)

        # get the averaged CTCF binding signal across all cancer cell lines for each group of CTCFs
        gained_signals =  signal_df.loc[gained_df.index.intersection(signal_df.index),cancercols].mean(axis=1)
        lost_signals =  signal_df.loc[lost_df.index.intersection(signal_df.index),cancercols].mean(axis=1)
        constitutive_signals =  signal_df.loc[constitutive_df.index.intersection(signal_df.index),cancercols].mean(axis=1)
        np.save(outdir+os.sep+'binding_level_on_cancer_datasets_{}.gained_signals'.format(cancertype),gained_signals)
        np.save(outdir+os.sep+'binding_level_on_cancer_datasets_{}.lost_signals'.format(cancertype),lost_signals)
        np.save(outdir+os.sep+'binding_level_on_cancer_datasets_{}.constitutive_signals'.format(cancertype),constitutive_signals)







if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
