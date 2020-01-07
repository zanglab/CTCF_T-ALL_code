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
# import return_CTCF_based_nearby_GeneInfo
import CTCF_TALL_modules_new

species = 'hg38'
gene_4k_bed = '/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'  
all_regions = association_with_regions.read_regions_from_bed_not_sep_strand(gene_4k_bed,species,expand = -1999)
region_start_list,region_end_list = association_with_regions.region_start_end_list(all_regions)


CTCF_domain_info_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo_GT100K_LT1M_EachSide.csv'
CTCF_domain_df = pd.read_csv(CTCF_domain_info_file,index_col=2)


def within_region_gene_ids(chr,region):
    # return ids of genes overlapped with given chr-region
    s = bisect.bisect_left(region_end_list[chr],region[0])
    e = bisect.bisect_right(region_start_list[chr],region[1])
    overlapped_IDs = []
    for i in np.arange(s,e):
        ID = all_regions[chr][i][-1].split('\t')[0]
        overlapped_IDs.append(ID)
    return overlapped_IDs  


def same_domain_overlapped_gene_ids(ctcf_id):
    # return ids of genes within same domain of CTCF 
    chr = CTCF_domain_df.loc[ctcf_id,'chr']
    domain_left = CTCF_domain_df.loc[ctcf_id,'domain_100k_1M_left']
    domain_right = CTCF_domain_df.loc[ctcf_id,'domain_100k_1M_right']
    #domain_len = CTCF_domain_df.loc[ctcf_id,'domain_len']
    if chr in all_regions:
    #if chr in chrs and (domain_len<5000000):
        overlapped_ids = within_region_gene_ids(chr,[domain_left,domain_right])
    else:
        overlapped_ids = []
    return overlapped_ids



def main():
    
    outdir = 'gene_CTCF_same_domain_GT100k_LT1M_prediction'
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


    for chr in chroms:
        #chr='chr21'
        print(chr)
        ctcf_ids = union_df[union_df['chr']==chr].index
        gene_ids = gene_domain_df[gene_domain_df['chr']==chr].index
        print('#CTCF:\t{}\t#genes:\t{}'.format(len(ctcf_ids),len(gene_ids)))
        combined_df = pd.DataFrame(index = ctcf_ids,columns = gene_ids)
        combined_df[:]=-1
        
        for ctcf_id in ctcf_ids:
            #ctcf_id = 537600
            same_domain_gene_ids = same_domain_overlapped_gene_ids(ctcf_id)#;print(same_domain_gene_ids)
            same_domain_gene_ids = set(same_domain_gene_ids).intersection(gene_ids)#;print(same_domain_gene_ids)
            combined_df.loc[ctcf_id][same_domain_gene_ids]=2        
            #print(combined_df.loc[537600,'U2AF1']);exit()
        combined_df.to_csv('{}/{}_{}.csv'.format(outdir,outdir,chr))
        #exit()


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
