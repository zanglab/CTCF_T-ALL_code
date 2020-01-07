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

def main():

    outdir='geneExpr_CTCFbinding_correlation'
    os.makedirs(outdir,exist_ok=True)
    
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
    'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

    # TSS position of each gene
    gene_domain_info_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f2_CTCF_binding_vs_GeneExpr_AllCancers/f2_ctcf_binding_GeneExpr_cor/f8_gene_domain_info/all_genes_domainInfo.csv'
    gene_domain_df = pd.read_csv(gene_domain_info_file,index_col=0)
    gene_domain_df = gene_domain_df.fillna(0)  
    # ucsc_df = association_with_genes.return_ucsc_df('hg38')
    gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
    gene_expr = pd.read_csv(gene_expr_file,index_col = 0)    
    
    # shared genes with TPM and TSS
    gene_collections = gene_domain_df.index.intersection(gene_expr.index)
    gene_domain_df = gene_domain_df.loc[gene_collections]
    # mid position of each CTCF binding site    
    union_df = CTCF_TALL_modules_new.return_occupancy_filtered()

    ctcf_gene_cor_dir = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f1_gene_CTCF_cor/f8_cor_ID_convert/CTCF_gene_cor_csv_ID_converted"
    
    for chr in chroms:
#         chr='chr4'
        print(chr)
        combined_df = pd.DataFrame()
        ctcf_ids = union_df[union_df['chr']==chr].index
        gene_ids = gene_domain_df[gene_domain_df['chr']==chr].index
        print('#CTCF:\t{}\t#genes:\t{}'.format(len(ctcf_ids),len(gene_ids)))
        for gene_id in gene_ids:
            #print(gene_id)
            ctcf_gene_info_file = ctcf_gene_cor_dir+os.sep+'{}_cor.csv'.format(gene_id)
            #if not os.path.isfile(ctcf_gene_info_file):
            #    print(chr,gene_id);#exit()
            ctcf_gene_info_df = pd.read_csv(ctcf_gene_info_file,header=None);
            ctcf_gene_info_df.index = ctcf_gene_info_df.index+1
            
            df_tmp = ctcf_gene_info_df.loc[ctcf_ids]
            df_tmp.columns = [gene_id]
            combined_df = pd.concat([combined_df,df_tmp],axis=1)
            #print(df_tmp);exit()
        combined_df = combined_df.round(4)
        #print(combined_df);exit()
        combined_df.to_csv('{}/geneExpr_CTCFbinding_correlation_{}.csv'.format(outdir,chr))
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
