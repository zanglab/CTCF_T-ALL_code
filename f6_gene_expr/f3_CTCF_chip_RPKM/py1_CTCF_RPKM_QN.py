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


def binding_signal_in_union():
    signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings.csv'
#     signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_head10000.csv'
    with open(signal_file) as signal_inf:
        signal_df = pd.read_csv(signal_inf,sep='\t',index_col=0)
    return signal_df



def main(indir):               

    outdir='f1_chip_RPKM_on_union_bindings'
    os.makedirs(outdir,exist_ok=True)
    
    matched_name_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f0_tables/RNA_ChIP_cor_match_IDs_201907.xlsx'
    matched_name_df = pd.read_excel(matched_name_file,index_col=0,sheet_name='ChIPseq_IDs') 
    print(matched_name_df.shape)

    signal_df = binding_signal_in_union()
    df_all = signal_df[matched_name_df['csv_names'].values]    
    df_all = df_all.fillna(0)
    df_all = df_all.round(2)
    df_all.to_csv('{}/RPKM_on_all_union_bindings_RPKM.csv'.format(outdir),sep='\t') 
    print(df_all.shape)

    ## ==== generate quantile normalized RPKM
    df_all = pd.read_csv('{}/RPKM_on_all_union_bindings_RPKM.csv'.format(outdir),sep='\t',index_col=0)
    df_new = pd.DataFrame()
    for matched_name_id in set(matched_name_df.index):
        csv_names = matched_name_df.loc[[matched_name_id]]['csv_names'].values
#         print(matched_name_id,csv_names)
#         print(matched_name_id,df_all[csv_names].shape)
        df_new = pd.concat([df_new,df_all[csv_names].mean(axis=1).rename(matched_name_id)],axis=1)
    
    
    df_new = df_new.round(2)
    df_new.to_csv('{}/RPKM_on_all_union_bindings_RPKM_CellType_combined.csv'.format(outdir)) 
        
    df = np.sqrt(df_new)
    df = quantile_normalization(df)
    df = df.round(2)
    df.to_csv('{}/RPKM_on_all_union_bindings_RPKM_CellType_combined_sqrt_QN.csv'.format(outdir)) 
        
    


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
