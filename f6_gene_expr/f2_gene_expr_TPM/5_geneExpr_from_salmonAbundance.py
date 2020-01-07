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
import association_with_genes
import association_with_regions
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
       


def main():
    
    outdir='f5_gene_expr_combined'
    os.makedirs(outdir,exist_ok=True)

    abundance_deseq2_file = 'f3_deseq_out/gene_level_abundance_TPM_all.csv'
    df = pd.read_csv(abundance_deseq2_file,index_col=0)
    #print(df);exit(0)    

    matched_name_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f0_tables/RNA_ChIP_cor_match_IDs_201907.xlsx'
    matched_name_df = pd.read_excel(matched_name_file,index_col=0,sheet_name='RNAseq_IDs') 
    print(matched_name_df.shape)

    df_new = pd.DataFrame()
    for matched_name_id in set(matched_name_df.index):
        csv_names = matched_name_df.loc[[matched_name_id]]['csv_names'].values
#         print(matched_name_id,csv_names)
#         print(matched_name_id,df_all[csv_names].shape)
        df_new = pd.concat([df_new,df[csv_names].mean(axis=1).rename(matched_name_id)],axis=1)

    df_new = df_new.round(2)
    df_new.to_csv('{}/geneExpr_replicates_combined_TPM.csv'.format(outdir))


    df_new = np.sqrt(df_new)
    df_new = df_new.round(2)
    df_new.to_csv('{}/geneExpr_replicates_combined_TPM_sqrt.csv'.format(outdir))



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
