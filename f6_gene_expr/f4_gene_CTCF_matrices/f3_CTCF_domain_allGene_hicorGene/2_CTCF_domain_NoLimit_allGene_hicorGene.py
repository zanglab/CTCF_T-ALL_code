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
import json,yaml


def return_dis_domain_cor_df(chrom):
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/fz_matrices'
    #chroms = return_chroms
    distance_file = pardir+os.sep+'gene_CTCF_pairwise_distances/gene_CTCF_pairwise_distances_{}.csv'.format(chrom)
    domain_info_file = pardir+os.sep+'gene_CTCF_same_domain_NoLimit_prediction/gene_CTCF_same_domain_NoLimit_prediction_{}.csv'.format(chrom)
    #domain_info_file = pardir+os.sep+'gene_CTCF_same_domain_GT100k_LT1M_prediction/gene_CTCF_same_domain_GT100k_LT1M_prediction_{}.csv'.format(chrom)
    cor_info_file = pardir+os.sep+'geneExpr_CTCFbinding_correlation/geneExpr_CTCFbinding_correlation_{}.csv'.format(chrom)
    with open(distance_file) as dis_inf,open(domain_info_file) as domain_inf, open(cor_info_file) as cor_inf:
        dis_df = pd.read_csv(dis_inf,index_col=0)
        domain_df = pd.read_csv(domain_inf,index_col=0)
        cor_df = pd.read_csv(cor_inf,index_col=0)
    dis_df = np.abs(dis_df)
    gene_collections = CTCF_TALL_modules_new.return_genome_gene_collections()
    dis_df = dis_df[dis_df.columns.intersection(gene_collections)]
    domain_df = domain_df[domain_df.columns.intersection(gene_collections)]
    cor_df = cor_df[cor_df.columns.intersection(gene_collections)]    
    return dis_df,domain_df,cor_df  


def main():
    
#     outdir = 'f1_CTCF_domain_GT100K_LT1M_allGene_hicorGene'
#     os.makedirs(outdir,exist_ok=True)

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
    'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
   
    ctcf_domain_gene_dic = {}
    for chr in chroms:
#         chr='chr21'
        print(chr)
        dis_df,domain_df,cor_df = return_dis_domain_cor_df(chr)
        promoter_df = cor_df[domain_df[(dis_df>=0)&(dis_df<=2000)]==2]
        promoter_df_hocor = promoter_df[np.power(promoter_df,2)>0.25]
        intra_domain_df = cor_df[domain_df==2]
        intra_domain_df_hocor = intra_domain_df[np.power(intra_domain_df,2)>0.25]
        for ctcf_id in dis_df.index:
            ctcf_domain_gene_dic[ctcf_id]={}
            ctcf_domain_gene_dic[ctcf_id]['promoter_genes']= promoter_df.loc[ctcf_id].dropna().index.to_list()
            ctcf_domain_gene_dic[ctcf_id]['promoter_genes_hicor']= promoter_df_hocor.loc[ctcf_id].dropna().index.to_list()
            ctcf_domain_gene_dic[ctcf_id]['intra_domain_genes']= intra_domain_df.loc[ctcf_id].dropna().index.to_list()
            ctcf_domain_gene_dic[ctcf_id]['intra_domain_genes_hicor']= intra_domain_df_hocor.loc[ctcf_id].dropna().index.to_list()

    with open('CTCF_domain_NoLimit_allGene_hicorGene.json','w') as json_file:
        json.dump(ctcf_domain_gene_dic, json_file)
#     with open('CTCF_domain_NoLimit_allGene_hicorGene.json', 'r') as json_file:
#         ctcf_domain_gene_dic = json.load(json_file)        

#     with open('CTCF_domain_NoLimit_allGene_hicorGene.yml','w') as yaml_outf:
#         yaml.dump(ctcf_domain_gene_dic, yaml_outf)
#     with open('CTCF_domain_NoLimit_allGene_hicorGene.yml', 'r') as yaml_inf:
#         ctcf_domain_gene_dic = yaml.load(yaml_inf)        





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
