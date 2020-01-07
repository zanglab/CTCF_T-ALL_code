'''
This file is used to get the TSS/domain/chr overlapped CTCFs for each gene
and return gene-CTCF cor given a geneID and a bounch of CTCFs
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new

def return_union_binding_bed_ID_convert():
#     union_CTCF='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f1_union_all_summits/f1_peak2000_datasets_union_summits/union_summits_fe4_width_150_sorted_natural.bed'
    union_CTCF='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f1_gene_CTCF_cor_old/f8_cor_ID_convert/union_summits_fe4_width_150_4thCol_NewNaturalID_LastCol_BedtoolSortID.bed'
    with open(union_CTCF) as union_inf:
        union_df = pd.read_csv(union_inf,sep='\t',index_col=6,header=None)
    union_df.columns = ['chr','start','end','nature_sort_ID','score','strand']
    return union_df



def main():

    outdir='CTCF_gene_cor_csv_ID_converted'
    os.makedirs(outdir,exist_ok=True)
    
    # TSS position of each gene
    gene_domain_info_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f2_CTCF_binding_vs_GeneExpr_AllCancers/f2_ctcf_binding_GeneExpr_cor/f8_gene_domain_info/all_genes_domainInfo.csv'
    gene_domain_df = pd.read_csv(gene_domain_info_file,index_col=0)
    gene_domain_df = gene_domain_df.fillna(0)  
    # ucsc_df = association_with_genes.return_ucsc_df('hg38')
    gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
    gene_expr = pd.read_csv(gene_expr_file,index_col = 0)    
    # shared genes with TPM and TSS
    gene_collections = gene_domain_df.index.intersection(gene_expr.index)
    
    # union df, with old bedtools sorted ID in index, new ID in col 'nature_sort_ID'    
    union_df_ID = return_union_binding_bed_ID_convert()
    
    ctcf_gene_cor_dir = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f1_gene_CTCF_cor_old/f7_cor_CTCF_binding_geneExpr/CTCF_gene_cor_csv"
    for gene_id in gene_collections[14000:16000]: 
        print(gene_id)
        ctcf_gene_info_file = ctcf_gene_cor_dir+os.sep+'{}_cor.csv'.format(gene_id)
        ctcf_gene_info_df = pd.read_csv(ctcf_gene_info_file,header=None);
        ctcf_gene_info_df.index = ctcf_gene_info_df.index+1
        df = pd.concat([union_df_ID,ctcf_gene_info_df],axis=1) 
        # sort by new ID, and save cor file 
        df = df[['nature_sort_ID',0]].sort_values(by=['nature_sort_ID'])[[0]]
        outfile = '{}/{}_cor.csv'.format(outdir,gene_id)
        with open(outfile,'w') as outf:
            for ii in df[0].values:
                outf.write('{:2f}\n'.format(ii))
        outf.close()
#         df = df.round(6)        
#         df.to_csv(outfile,index=False,header=None)
#         exit()




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--geneid', action = 'store', type = str,dest = 'geneid', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
