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

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')

def sqrt_divid_med_compr(c1,c2):
    med1 = np.median(c1)
    med2 = np.median(c2);print(med1,med2)
    return np.sqrt(c1/med1)-np.sqrt(c2/med2)

def quantile_normalization(dataframe):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = dataframe
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized

def log_quantile_normalization(df):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = np.log2(df+1)
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized
             
def col_row_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    df = df.subtract(df.mean(axis=1),axis=0)
    return df


def main():

   # if species in species_chroms.keys():
   #     chroms = species_chroms[species]
    indir = 'binding_csv'
    outdir = 'binding_changes_QN'
    os.makedirs(outdir,exist_ok=True)
    
    
    suffixs = ['e200','e500','binding_site']
    treats = ['dmso','gsi_3d','gsi_wo']
    
    for suffix in suffixs:
        # for each binding region, get the binding change score
        df = pd.DataFrame()
        for treat in treats:
            file = indir+os.sep+'{}_ATAC_RPKM_{}.csv'.format(treat,suffix)
            with open(file) as inf:
                treat_df = pd.read_csv(inf,sep='\t',header=None,index_col=0)
                treat_df.columns = ['{}'.format(treat)]
            df = pd.concat([df,treat_df],axis=1)

        print(suffix)  
        df1 = log_quantile_normalization(df)     
        df1['{}_vs_{}'.format(treats[1],treats[0])] = df1[treats[1]]-df1[treats[0]]
        df1['{}_vs_{}'.format(treats[2],treats[1])] = df1[treats[2]]-df1[treats[1]]
        df1['{}_vs_{}'.format(treats[2],treats[0])] = df1[treats[2]]-df1[treats[0]]
        
        df1.to_csv(outdir+os.sep+'ATAC_RPKM_changes_on_Union_{}.csv'.format(suffix))
        

    #### try to use only ATAC peak overlapped regions
        
        peak_overlapped_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f11_GSI_shCTCF/f5_GSI_ATAC_chip_20191029/f2_results_from_fq_mapping/f1_accessibility_signal_csv/union_binding_mid_GSI_ATAC_peak_overlapped.bed'
        peak_overlapped_df =  pd.read_csv(peak_overlapped_file,sep='\t',index_col=3,header=None)  
        df = df.loc[peak_overlapped_df.index].dropna()
        df1 = log_quantile_normalization(df)     
        print('overlapped',suffix)        
        df1['{}_vs_{}'.format(treats[1],treats[0])] = df1[treats[1]]-df1[treats[0]]
        df1['{}_vs_{}'.format(treats[2],treats[1])] = df1[treats[2]]-df1[treats[1]]
        df1['{}_vs_{}'.format(treats[2],treats[0])] = df1[treats[2]]-df1[treats[0]]
        
        df1.to_csv(outdir+os.sep+'ATAC_RPKM_changes_on_Union_peak_overlapped_union_{}.csv'.format(suffix))


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
