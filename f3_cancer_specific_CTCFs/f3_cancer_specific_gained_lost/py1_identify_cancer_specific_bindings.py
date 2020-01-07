import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new


def identification_of_gained_bindings(df,cancertype):
    
    # ==== set threshold for different cancer types ====
    all_total,normal_total,cancer_total = df.iloc[0,:]['all_total'], df.iloc[0,:]['normal_total'], df.iloc[0,:]['cancer_total']
    all_thre,normal_thre,cancer_thre,all_pthre,normal_pthre = int(all_total*0.2), int(normal_total*0.0), max(int(cancer_total*0.5),2), 0.01, .05
    
    print('\n{}\t gained\t all-total:{}\t normal-total:{}\t cancer-total:{}'.format(cancertype,all_total,normal_total,cancer_total))
    print('{}\t gained\t all-thre:{}\t normal-thre:{}\t cancer-thre:{}\t pthre:{}\n'.format(cancertype,all_thre,normal_thre,cancer_thre,all_pthre))
    
    df = df[df['all_peaks'] <= all_thre]
    df = df[df['normal_peaks'] <= normal_thre]
    df = df[df['cancer_peaks'] >= cancer_thre]
    
    df = df[df['cancer_vs_other_stats']>0]
    df = df[df['cancer_vs_other_pvalue']<=all_pthre]
    df = df[df['cancer_vs_other_adjp']<=all_pthre]
    
    df = df[df['cancer_vs_normal_stats']>0]
#     df = df[df['cancer_vs_normal_pvalue']<=normal_pthre]
#     df = df[df['cancer_vs_normal_adjp']<=normal_pthre]
    
    return df
    


def identification_of_lost_bindings(df,cancertype):

    all_total,normal_total,cancer_total = df.iloc[0,:]['all_total'], df.iloc[0,:]['normal_total'], df.iloc[0,:]['cancer_total']
    all_thre,normal_thre,cancer_thre,all_pthre,normal_pthre = int(all_total*0.7), max(2,int(normal_total*0.5)), int(cancer_total*0.2), .1, .1
    
    print('{}\t lost\t all-total:{}\t normal-total:{}\t cancer-total:{}'.format(cancertype,all_total,normal_total,cancer_total))
    print('{}\t lost\t all-thre:{}\t normal-thre:{}\t cancer-thre:{}\t pthre:{}\n'.format(cancertype,all_thre,normal_thre,cancer_thre,all_pthre))

    df = df[df['all_peaks'] >= all_thre]
    df = df[df['normal_peaks'] >= normal_thre]
    df = df[df['cancer_peaks'] <= cancer_thre]
    
    df = df[df['cancer_vs_other_stats']<0]
#     df = df[df['cancer_vs_other_pvalue']<=all_pthre]
#     df = df[df['cancer_vs_other_adjp']<=all_pthre]
    
    df = df[df['cancer_vs_normal_stats']<0]
#     df = df[df['cancer_vs_normal_pvalue']<=normal_pthre]
#     df = df[df['cancer_vs_normal_adjp']<=normal_pthre]
    
    return df

    


def main():

    outdir = 'f1_cancer_specific_binding'
    os.makedirs(outdir,exist_ok=True)
    
    # keep only bindings with occupancy score ≥3
    occupancy_filtered_df = CTCF_TALL_modules_new.return_occupancy_filtered()
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
#     cancertypes=['AML']
    for cancertype in cancertypes: 
        all_feature_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f2_combined_sep_cancertype/{}_binding_features.csv'.format(cancertype)
        with open(all_feature_file) as all_feature_inf:
            df = pd.read_csv(all_feature_inf,index_col=0)
        
        # ==== gained bindings for each cancer type ====
        gained_bindings = identification_of_gained_bindings(df,cancertype)
        gained_df = df.loc[gained_bindings.index]
        
        # ==== lost bindings for each cancer type ====
        lost_bindings = identification_of_lost_bindings(df,cancertype)
        lost_df = df.loc[lost_bindings.index]
      
        # ==== save gained/lost bindings ====
        gained_df.to_csv(outdir+os.sep+'{}_gained.csv'.format(cancertype))
        lost_df.to_csv(outdir+os.sep+'{}_lost.csv'.format(cancertype))
    
        gained_df = gained_df[['chr','start','end','#motifs']];gained_df['strand']='+'
        gained_df.insert(3,'id',gained_df.index)
        gained_df.to_csv(outdir+os.sep+'{}_gained.bed'.format(cancertype),sep='\t',header=None,index=False)
    
        lost_df = lost_df[['chr','start','end','#motifs']];lost_df['strand']='+'
        lost_df.insert(3,'id',lost_df.index)
        lost_df.to_csv(outdir+os.sep+'{}_lost.bed'.format(cancertype),sep='\t',header=None,index=False)
    
        print(cancertype);print('#gained',gained_df.shape);print('#lost',lost_df.shape,'\n====')
#         exit()
    




   
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

