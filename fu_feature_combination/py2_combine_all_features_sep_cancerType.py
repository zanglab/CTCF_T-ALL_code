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
import CTCF_TALL_modules_new


chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

def add_patient_peaks(df):

    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f9_TALL_patient_peakss/f2_patient_overlap_counting/union_binding_patient_peak_info.csv'  ### need to be changed
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=3) ### need to be changed
    #print(add_df,add_df.columns);exit()
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['PD31_peak', 'PD9_peak', 'all_patients_peaks'] # add_df.columns ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df


def add_dynamic_notch_intra_domain_co_occurrence(df):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f1_ctcf_notch_enrichment/f1_notch_myc_enrichment/dynamic_notch_domain_co_occurrence.csv'
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=4,sep='\t',low_memory=False) ### need to be changed
    add_df.columns=['#chrom','start','end','if_intra_domain_dynamic_notch','score','strand','dynamic_notch_center']
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['if_intra_domain_dynamic_notch','dynamic_notch_center'] # add_df.columns ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def add_cutll1_notch_intra_domain_co_occurrence_new(df):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f5_CTCF_notch_targetGene/f1_notch_coocurrence_with_center_pos/cutll1_w4h_notch_domain_co_occurrence.csv'
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=4,sep='\t',low_memory=False) ### need to be changed
    add_df.columns=['#chrom','start','end','if_intra_domain_notch','score','strand','notch_center']
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['if_intra_domain_notch','notch_center'] # add_df.columns ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def add_cutll1_notch_intra_domain_co_occurrence(df):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f1_ctcf_notch_enrichment/f1_notch_myc_enrichment/cutll1_notch_domain_co_occurrence.csv'
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=3,sep='\t',header=None) ### need to be changed
    add_df.columns=['#chrom','start','end','score','strand','if_intra_domain_notch']
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['if_intra_domain_notch'] # add_df.columns ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df
        
def add_cutll1_MYC_intra_domain_co_occurrence(df):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f10_TF_binding/f3_TALL_Notch/f1_ctcf_notch_enrichment/f1_notch_myc_enrichment/cutll1_MYC_domain_co_occurrence.csv'
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=3,sep='\t',header=None) ### need to be changed
    add_df.columns=['#chrom','start','end','score','strand','if_intra_domain_MYC']
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['if_intra_domain_MYC'] # add_df.columns ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def add_DNA_methylation_cov5(df,cancertype):
    cancer_type_matched = {'LUAD':'A549_cov5_vs_lung_cov5_WGBS','CRC':'HCT116_cov5_vs_sigmoid_colon_cov5','BRCA':'MCF7_cov5_vs_breast_cov5',\
                           'T-ALL':'Jurkat_vs_A6010','PD9':'PD9_cov5_vs_A6010','PD30':'PD31_cov5_vs_A6010'}
    if cancertype in cancer_type_matched.keys():
        add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f2_union_binding_methylation_compr/f2_union_binding_differential_DNA_methylation/{}.csv'.format(cancer_type_matched[cancertype])
        with open(add_file) as add_inf:
            add_df = pd.read_csv(add_inf,index_col=0) ### need to be changed
        assert len(df.index.difference(add_df.index))==0
        add_df = add_df.loc[df.index]
        add_columns = add_df.columns ### need to be changed
        df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def add_deg_hicor_info(df,cancertype):
    if cancertype=='T-ALL':
        cancertype='T_ALL'
    if cancertype=='PRAD_TissueAdded':
        cancertype='PRAD'
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f1_ctcf_domain_GT100k_LT1M_DEG/{}_domain_GT100k_LT1M_DEG_log2FC0.585_padj1e-3.csv'.format(cancertype)  ### need to be changed
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f1_ctcf_domain_GT100k_LT1M_DEG/{}_domain_GT100k_LT1M_DEG_log2FC1_padj1e-5.csv'.format(cancertype)  ### need to be changed
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f2_ctcf_domain_NoLimit_DEG/{}_domain_NoLimit_DEG_log2FC0.585_padj1e-3.csv'.format(cancertype)  ### need to be changed
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f7_panCancer_deg_figs/f2_ctcf_domain_NoLimit_DEG/{}_domain_NoLimit_DEG_log2FC1_padj1e-5.csv'.format(cancertype)  ### need to be changed
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=0,low_memory=False) ### need to be changed
    #print(add_df.columns);exit()
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = add_df.columns ### need to be changed
    '''
    add_df.columns = [
       promoter_allgenes,promoter_allgenes_pos,promoter_allgenes_hicor,promoter_allgenes_hicor_pos,
       promoter_upgenes,promoter_upgenes_pos,promoter_upgenes_hicor,promoter_upgenes_hicor_pos,promoter_upgenes_names,promoter_upgenes_hicor_names,
       promoter_dngenes,promoter_dngenes_pos,promoter_dngenes_hicor,promoter_dngenes_hicor_pos,promoter_dngenes_names,promoter_dngenes_hicor_names,
       domain_allgenes,domain_allgenes_pos,domain_allgenes_hicor,domain_allgenes_hicor_pos,
       domain_upgenes,domain_upgenes_pos,domain_upgenes_hicor,domain_upgenes_hicor_pos,domain_upgenes_names,domain_upgenes_hicor_names,
       domain_dngenes,domain_dngenes_pos,domain_dngenes_hicor,domain_dngenes_hicor_pos,domain_dngenes_names,domain_dngenes_hicor_names ]
    '''
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def add_binding_signal(df,cancertype):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f2_cancerType_binding_signal/f1_cancertype_signal_ttest_fdr/{}_signals_ttest_fdr.csv'.format(cancertype)  ### need to be changed
    if cancertype in ['PD9','PD30']:
        add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f2_cancerType_binding_signal/f1_cancertype_signal_ttest_fdr/T-ALL_signals_ttest_fdr.csv'.format(cancertype)  ### need to be changed
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,sep='\t',index_col=0) ### need to be changed
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['cancer_mean','normal_mean','ctrl_mean','cancer_vs_other_stats','cancer_vs_other_pvalue','cancer_vs_other_adjp',\
    'normal_vs_other_stats','normal_vs_other_pvalue','normal_vs_other_adjp',\
    'cancer_vs_normal_stats','cancer_vs_normal_pvalue','cancer_vs_normal_adjp'] ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)
    return df

def add_binding_occupancy(df,cancertype):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f1_cancerType_occupancy/f1_cancerType_binding_occurrence/{}_binding_occurrence.csv'.format(cancertype)  ### need to be changed
    if cancertype in ['PD9','PD30']:
        add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f1_cancerType_occupancy/f1_cancerType_binding_occurrence/T-ALL_binding_occurrence.csv'.format(cancertype)  ### need to be changed
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=0,sep=',') ### need to be changed
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['cancer_peaks','normal_peaks','all_peaks','cancer_total','normal_total','all_total'] ### need to be changed
    df = pd.concat([df,add_df[add_columns]],axis=1)
    return df

def add_domain_position(df,cancertype):
    add_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo_GT100K_LT1M_EachSide.csv'
    with open(add_file) as add_inf:
        add_df = pd.read_csv(add_inf,index_col=2)
    assert len(df.index.difference(add_df.index))==0
    add_df = add_df.loc[df.index]
    add_columns = ['domain_100k_1M_left','domain_100k_1M_right']
    df = pd.concat([df,add_df[add_columns]],axis=1)#;print(df);exit()
    return df

def combine_sep_cancer_types(cancertype):
    outdir = 'f2_combined_sep_cancertype'
    os.makedirs(outdir,exist_ok=True)
    
    df = CTCF_TALL_modules_new.return_occupancy_filtered()
    df = add_domain_position(df,cancertype)
    
    # == for 'PD9','PD30', use same data as T-ALL
    df = add_binding_occupancy(df,cancertype) # cancer total: 7, normal total: 3, all total: 584
    df = add_binding_signal(df,cancertype)

    # == matched data for each cancer type
    df = add_deg_hicor_info(df,cancertype)
    df = add_DNA_methylation_cov5(df,cancertype)

    if cancertype=='T-ALL':
#         df = add_cutll1_notch_intra_domain_co_occurrence(df)
        df = add_cutll1_notch_intra_domain_co_occurrence_new(df)
        df = add_dynamic_notch_intra_domain_co_occurrence(df)
        df = add_cutll1_MYC_intra_domain_co_occurrence(df)
        df = add_patient_peaks(df)
    df.to_csv(outdir+os.sep+'{}_binding_features.csv'.format(cancertype))
    exit()

def main():

    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded','PD9','PD30']
#     cancertypes=['PD9','PD30']
    for cancertype in cancertypes:
        print(cancertype)
        combine_sep_cancer_types(cancertype)
        exit()


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

