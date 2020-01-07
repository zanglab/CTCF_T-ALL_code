import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
from stats_pvalues import irwin_hall_cdf

cancertype_myname = {'BRCA':'BRCA','COAD':'CRC','LUAD':'LUAD'}


def irwin_hall_return_pvalue(df):
    # for each sample/col, get the rank (smaller is better) of each item (row)
    # calculate the rank sum of all items across all samples
    df_rank_sum = df.rank(ascending=False).sum(axis=1)/df.shape[0]
    df['rank_sum'] = df_rank_sum
#     df['irwin_hall_pvalue']=[irwin_hall_cdf(i,df.shape[1]) for i in df_rank_sum]
    return df
    
    
def main(infile):
    
    outdir = 'f2_ATAC_sample_irwin_hall'
    os.makedirs(outdir,exist_ok=True)
    
    
    # selected cancer-specific ATAC info (overlapped with cancer specific CTCFs)
    atac_sig_file_dir = 'f1_ATAC_sig_caseID_differential_score'
    cancertypes = ['BRCA','CRC','LUAD','PRAD','PRAD_TissueAdded']
    cancertype_match_names={'BRCA':'BRCA','CRC':'COAD','LUAD':'LUAD','PRAD':'PRAD','PRAD_TissueAdded':'PRAD'}
    
    for name_type in cancertypes:
        for binding_type in ['gained','lost']:
            sig_file = atac_sig_file_dir+os.sep+'{}_{}_caseID_ATAC_sig_diff_score.csv'.format(name_type,binding_type)
            #print(sig_file)
            with open(sig_file) as sig_inf:
                sig_df = pd.read_csv(sig_inf,sep='\t')
#             print(sig_df);exit()
            sig_df = sig_df.iloc[:,7:]
            sig_df = np.transpose(sig_df)#;print(sig_df)
            sig_df_p = irwin_hall_return_pvalue(sig_df) # rank patient samples for each region/col
            sig_df_p = sig_df_p.sort_values(by=['rank_sum'],ascending=True)
            sig_df_p.to_csv(outdir+os.sep+'{}_{}_irwin_hall.csv'.format(name_type,binding_type))
            
            cancer_type_index = [i for i in sig_df.index if re.search(cancertype_match_names[name_type],i)]
            sig_df_p = irwin_hall_return_pvalue(sig_df.loc[cancer_type_index])
            sig_df_p = sig_df_p.sort_values(by=['rank_sum'],ascending=True)
            sig_df_p.to_csv(outdir+os.sep+'{}_{}_irwin_hall_cancertype_ID.csv'.format(name_type,binding_type))

                
             


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
