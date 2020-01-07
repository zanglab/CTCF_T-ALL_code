'''
for each CTCF binding
    for TSS group, intra_domain gene group, return
    
        # all_genes, hicor_all_genes, hicor_pos_all_genes
        # up_genes, hicor_up_genes, hicor_pos_up_genes
        # dn_genes, hicor_dn_genes, hicor_pos_dn_genes

'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#import return_deg
import CTCF_TALL_modules_new


def return_dis_domain_cor_df(chrom):
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/fz_matrices'
    #chroms = return_chroms
    distance_file = pardir+os.sep+'gene_CTCF_pairwise_distances/gene_CTCF_pairwise_distances_{}.csv'.format(chrom)
    #domain_info_file = pardir+os.sep+'gene_CTCF_same_domain_NoLimit_prediction/gene_CTCF_same_domain_NoLimit_prediction_{}.csv'.format(chrom)
    domain_info_file = pardir+os.sep+'gene_CTCF_same_domain_GT100k_LT1M_prediction/gene_CTCF_same_domain_GT100k_LT1M_prediction_{}.csv'.format(chrom)
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
    
    

def abstract_hicor_info(df,return_df,flag,cor_thre=0.5,hicor_gene_name_saved=False):      

    return_df['{}'.format(flag)] = np.abs(np.sign(df)).sum(axis=1)
    return_df['{}_pos'.format(flag)] = np.sign(df).clip(0,1).sum(axis=1)
    #print(df);exit()
    hicor_df = df[(df>cor_thre)|(df<-1*cor_thre)]  
    return_df['{}_hicor'.format(flag)] = np.abs(np.sign(hicor_df)).sum(axis=1)
    return_df['{}_hicor_pos'.format(flag)] = np.sign(hicor_df).clip(0,1).sum(axis=1)
    
    if hicor_gene_name_saved:
        return_df['{}_names'.format(flag)] = df.apply(lambda x: ','.join(x[x.notnull()].index),axis=1)
        return_df['{}_hicor_names'.format(flag)] = hicor_df.apply(lambda x: ','.join(x[x.notnull()].index),axis=1)
              
    return return_df


def return_ctcf_gene_cor_info(dataType,outdir,up_genes,dn_genes,flag):

    info_df = pd.DataFrame()
    chroms = CTCF_TALL_modules_new.return_chroms()
#     chroms=['chr21','chr22']
    for chrom in chroms:
        dis_df,domain_df,cor_df = return_dis_domain_cor_df(chrom)
        chr_info_df = pd.DataFrame()
        
        # look for if promoter genes are diffExpr, and if have hicor with CTCF
        promoter_cor_df = cor_df[domain_df[(dis_df>=0)&(dis_df<=2000)]==2]
        chr_info_df = abstract_hicor_info(promoter_cor_df,chr_info_df,'promoter_allgenes')
        # promoter up-regulated genes
        promoter_cor_df_upgenes = promoter_cor_df[promoter_cor_df.columns.intersection(up_genes)]
        chr_info_df = abstract_hicor_info(promoter_cor_df_upgenes,chr_info_df,'promoter_upgenes',hicor_gene_name_saved=True)
        # promoter down-regulated genes
        promoter_cor_df_dngenes = promoter_cor_df[promoter_cor_df.columns.intersection(dn_genes)]
        chr_info_df = abstract_hicor_info(promoter_cor_df_dngenes,chr_info_df,'promoter_dngenes',hicor_gene_name_saved=True)

        # look for if intra-domain genes are diffExpr, and if have hicor with CTCF
        domain_cor_df = cor_df[domain_df==2]
        chr_info_df = abstract_hicor_info(domain_cor_df,chr_info_df,'domain_allgenes')
        # intra-domain up-regulated genes
        domain_cor_df_upgenes = domain_cor_df[domain_cor_df.columns.intersection(up_genes)]
        chr_info_df = abstract_hicor_info(domain_cor_df_upgenes,chr_info_df,'domain_upgenes',hicor_gene_name_saved=True)
        # intra-domain down-regulated genes
        domain_cor_df_dngenes = domain_cor_df[domain_cor_df.columns.intersection(dn_genes)]
        chr_info_df = abstract_hicor_info(domain_cor_df_dngenes,chr_info_df,'domain_dngenes',hicor_gene_name_saved=True)
        
        info_df = pd.concat([info_df,chr_info_df])
    info_df.to_csv(outdir+os.sep+'{}_domain_GT100k_LT1M_DEG_{}.csv'.format(dataType,flag))


def return_cancertype_deg_TCGA(dataType,fc_thre=0.585,padj_thre=0.01):
    # deg identified by final collection, use for ctcf-gene correlation
    dataTypes = ['BRCA','COAD','LUAD','PRAD']
    
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f8_panCancer_DEG_TCGA/tcga_DESeq_Results'
    deseq2_file = pardir+os.sep+'/{}/DESeq2ResultsGeneSymbols.csv'.format(dataType)
    if dataType=='CRC':
        deseq2_file = pardir+os.sep+'/{}/DESeq2ResultsGeneSymbols.csv'.format(dataType)
    
    with open(deseq2_file) as deseq2_inf:
        deseq2_df = pd.read_csv(deseq2_inf,index_col=8) 
    deseq2_df = deseq2_df[deseq2_df.index.notnull()]
#     print(deseq2_df);exit()
    
    deseq2_up_genes = deseq2_df[(deseq2_df['log2FoldChange']>fc_thre)&(deseq2_df['padj']<padj_thre)].index.values
    deseq2_down_genes = deseq2_df[(deseq2_df['log2FoldChange']<-1*fc_thre)&(deseq2_df['padj']<padj_thre)].index.values 
    return deseq2_up_genes,deseq2_down_genes


def main():

    outdir = 'f1_ctcf_domain_GT100k_LT1M_DEG'
    os.makedirs(outdir,exist_ok=True)
        
    # ==== datatype for DEG
    dataTypes = ['AML','BRCA','CRC','LUAD','PRAD','T_ALL','T_ALL_CUTLL1','T_ALL_Jurkat','PD9','PD30']
    dataTypes = ['BRCA','CRC','LUAD','PRAD']
    
    fc_thre=0.585
    padj_thre=0.001
    print('Threshold:\t',fc_thre,padj_thre)
    for dataType in dataTypes:
        up_genes,dn_genes = return_cancertype_deg_TCGA(dataType,fc_thre,padj_thre)
        print(dataType,len(up_genes),len(dn_genes))
        return_ctcf_gene_cor_info(dataType,outdir,up_genes,dn_genes,'log2FC0.585_padj1e-3')
#         exit()
        

    fc_thre=1
    padj_thre=0.00001
    print('Threshold:\t',fc_thre,padj_thre)
    for dataType in dataTypes:
        up_genes,dn_genes = return_cancertype_deg_TCGA(dataType,fc_thre,padj_thre)
        print(dataType,len(up_genes),len(dn_genes))
        return_ctcf_gene_cor_info(dataType,outdir,up_genes,dn_genes,'log2FC1_padj1e-5')



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
