import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
#import association_with_regions
from get_reads_positions import reads_positions
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"

def quantile_normalization(dataframe):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = dataframe
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized


   
def main():

    outdir = 'f2_combined_RPKM_csv'
    os.makedirs(outdir,exist_ok=True)

    df = pd.DataFrame()
    
    for gsm in ['GSM1155964','GSM1155965','GSM1155966','GSM1155967','GSM1155968','GSM1155969','GSM2411156','GSM2411157','GSM2411158']:
        csv_file = 'f1_TALL_ATAC_RPKM_on_union/{}.txt'.format(gsm)
        with open(csv_file) as csv_inf:
            sample_df = pd.read_csv(csv_inf,sep="\t",index_col=0,header=None)
            sample_df.columns=[gsm]#;print(sample_df);exit()
        df = pd.concat([df,sample_df],axis=1)
    df = df.round(2)     
    df.to_csv(outdir+os.sep+'Jurkat_CD4_ATAC_combined_RPKM_raw.csv')
    
    df = np.log2(df+1)
    df = quantile_normalization(df) 
    df = df.round(2)     
    df.to_csv(outdir+os.sep+'Jurkat_CD4_ATAC_combined_RPKM_log2_QN.csv')

    



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
