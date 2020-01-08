import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=.7)
#sns.set_style("whitegrid", {'axes.grid' : False})
#import re,bisect
import json

####
# basic command lines
####

def basic_command_lines():
    
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()

    # intra-domain genes for each CTCF
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()

    # ==== genome wide genes    
    all_genes = CTCF_TALL_modules_new.return_genome_gene_collections()
    
    fcthre=np.log2(1.2)
    pthre=0.05    
    # ==== CUTLL1 up-genes
    cutll1_all,cutll1_upgenes,cutll1_dngenes = CTCF_TALL_modules_new.return_CUTLL1_vs_CD4_deg(fcthre,pthre)  
    tall_all,tall_upgenes,tall_dngenes = CTCF_TALL_modules_new.return_TALL_vs_CD4_deg(fcthre,pthre)
    # ==== NOTCH1 targets
    notch_targets = CTCF_TALL_modules_new.return_Notch1_target_genes()    
    # ==== shCTCF down regulated genes
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = CTCF_TALL_modules_new.return_CUTLL1_shCTCF_deg(fcthre,pthre)    
    # ==== GSI down genes
    gsi_all,gsi_upgenes,gsi_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_deg(fcthre,pthre)  
    # ==== GSI washout up genes
    gsi_wo_all,gsi_wo_upgenes,gsi_wo_dngenes = CTCF_TALL_modules_new.return_CUTLL1_GSI_wo_deg(fcthre,pthre)  


##
def return_ctcf_intra_domain_genes(binding_type,colname):
    # intra domain genes of CTCF bindings
    ctcf_domain_gene_dic = CTCF_TALL_modules_new.return_union_domain_gene()
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    
    tall_gained_intra_domain_genes = np.array([])
    for gained_id in gained_df.index:
        domain_genes = ctcf_domain_gene_dic[str(gained_id)]['intra_domain_genes']
        tall_gained_intra_domain_genes = np.append(tall_gained_intra_domain_genes,domain_genes)
    return tall_gained_intra_domain_genes


######################################
### basic collection/data
######################################

def return_chroms():
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    return chroms

def return_genome_gene_collections():
    # TSS position of each gene
    gene_domain_info_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f2_CTCF_binding_vs_GeneExpr_AllCancers/f2_ctcf_binding_GeneExpr_cor/f8_gene_domain_info/all_genes_domainInfo.csv'
    gene_domain_df = pd.read_csv(gene_domain_info_file,index_col=0)
    gene_domain_df = gene_domain_df.fillna(0)  
    # ucsc_df = association_with_genes.return_ucsc_df('hg38')
    gene_expr_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f2_gene_expr_TPM/f5_gene_expr_combined/geneExpr_replicates_combined_TPM_sqrt.csv"
    gene_expr = pd.read_csv(gene_expr_file,index_col = 0)    
    # genes with TSS 4kb info
    gene_4k_bed = '/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'  
    gene_4k_df = pd.read_csv(gene_4k_bed,sep='\t',index_col=3,header=None)
    
    # shared genes with TPM and TSS
    gene_collections = gene_domain_df.index.intersection(gene_expr.index).intersection(gene_4k_df.index)
    return gene_collections
    

def return_collection_df():
    ctcf_collection_file='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f0_data_collection_new/GEO_collection_tables/CTCF_peak2000_20190704.xlsx'
    collection_df = pd.read_excel(ctcf_collection_file,index_col=0)
    return collection_df

def peak_overlap_in_union():
    peak_overlap_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f1_each_data_peak_overlap_union/union_summits_EachDataOverlapInfo.csv'
#     peak_overlap_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f1_each_data_peak_overlap_union/union_summits_EachDataOverlapInfo_head10000.csv'
    with open(peak_overlap_file) as peak_overlap_inf:
        peak_overlap_df = pd.read_csv(peak_overlap_inf,sep="\t",index_col = 0)  
    return peak_overlap_df

def binding_signal_in_union():
    signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized.csv'
#     signal_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f3_signals_on_union_CTCFs/f2_signals_on_union_bindings/signals_RPKM_on_all_CTCF_bindings_QuantileNormalized_head10000.csv'
    with open(signal_file) as signal_inf:
        signal_df = pd.read_csv(signal_inf,sep='\t',index_col=0)
    return signal_df


######################################
### Union bindings
######################################

def return_union_binding_bed():
    union_CTCF='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f1_union_all_summits/f1_peak2000_datasets_union_summits/union_summits_fe4_width_150_sorted_natural.bed'
    with open(union_CTCF) as union_inf:
        union_df = pd.read_csv(union_inf,sep='\t',index_col=3,header=None)
    union_df.columns = ['chr','start','end','score','strand']
    return union_df

def return_union_binding_df():
    union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding.csv"
    with open(union_file) as union_inf:
        union_df = pd.read_csv(union_inf,sep=',',index_col=3)
    return union_df

def return_occupancy_filtered():
    union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.csv"
    with open(union_file) as union_inf:
        union_df = pd.read_csv(union_inf,sep=',',index_col=3)
    return union_df

def return_union_domain_gene():
    # ==== domain gene for each CTCF
    json_file_dir= '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f3_CTCF_domain_allGene_hicorGene'
#     with open(json_file_dir+os.sep+'CTCF_domain_NoLimit_allGene_hicorGene.json', 'r') as json_file:
    with open(json_file_dir+os.sep+'CTCF_domain_GT100K_LT1M_allGene_hicorGene.json', 'r') as json_file:
        ctcf_domain_gene_dic = json.load(json_file)        
    return ctcf_domain_gene_dic
    
    

def return_CTCF_domain_boundary():
    CTCF_domain_info_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo_GT100K_LT1M_EachSide.csv'
    with open(CTCF_domain_info_file) as CTCF_domain_info_inf:
        CTCF_domain_df = pd.read_csv(CTCF_domain_info_inf,index_col=2)
    CTCF_domain_df = CTCF_domain_df.fillna(0)
    return CTCF_domain_df


######################################
### Union bindings tmp
######################################

def return_union_binding_with_motif():
    motif_bed = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f4_motif_on_union_CTCFs/f1_union_binding_motif/union_CTCF_with_motif.bed'  
    with open(motif_bed) as motif_inf:
        motif_df = pd.read_csv(motif_inf,sep='\t',index_col=3)
    return motif_df

def return_occupancy_df():
    input_bed = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f2_occupancy_score/union_CTCF_occupancy_score.csv'  
    with open(input_bed) as input_inf:
        input_df = pd.read_csv(input_inf,index_col=3,sep='\t')
    return input_df


######################################
### Cancer specific bindings/data
######################################

def cancer_specific_cancer_normal_gsmID(cancertype):
    # cancer and normal GSM IDs for each cancer type
    cancerType_data = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f0_infile/CancerTypes_datasets_201907.xlsx'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    cancertype_df = pd.read_excel(cancerType_data,index_col=0,sheetname=cancertype) 
    cancercols = cancertype_df['CancerData'].dropna().values
    normalcols = cancertype_df['NormalData'].dropna().values
    return cancercols,normalcols

def return_cancer_specific_binding(cancertype):
    # == cancer specific gained/lost bindings ====
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f3_cancer_specific_gained_lost/f1_cancer_specific_binding'
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f3_cancer_specific_binding'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded','PD9','PD30']
    gained_file = pardir+os.sep+'{}_gained.csv'.format(cancertype)
    lost_file = pardir+os.sep+'{}_lost.csv'.format(cancertype)
    with open(gained_file) as gained_inf, open(lost_file) as lost_inf:
        gained_df = pd.read_csv(gained_inf,index_col=0)
        lost_df = pd.read_csv(lost_inf,index_col=0)
    return gained_df,lost_df

def return_cancer_specific_combined_features(cancertype):
    all_feature_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/fu_feature_combination/f2_combined_sep_cancertype/{}_binding_features.csv'.format(cancertype)
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded','PD9','PD30']
    with open(all_feature_file) as all_feature_inf:
        combined_df = pd.read_csv(all_feature_inf,index_col=0,low_memory=False)
    return combined_df


######################################
### Constitutive bindings
######################################

def return_constitutive_df():
    constitutive_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.csv'
    with open(constitutive_file) as constitutive_inf:
        constitutive_df = pd.read_csv(constitutive_inf,sep=',',index_col = 3)
    return constitutive_df

# def return_const_selected_ctrl():
#     filter_const_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f3_dynamic_HiC/f2_Hic_Pro_dynamic_compr/f1_compr_normalized_view_region_signals/f7_proper_ctrl_constitutive_selection/proper_const_ctrl.csv'
#     filter_const_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f2_bindings_signals_fdr/f3_proper_ctrl_constitutive_selection/const_NONoverlapped_panCancer_fdr_LT_001_TALL_log2fcLT1.bed'
#     with open(filter_const_file) as filter_const_inf:
#         filter_const_df = pd.read_csv(filter_const_inf,sep='\t',index_col=3,header=None)
#     return filter_const_df.index


# def return_cancertype_const_selected_ctrl(cancertype):
#     filter_const_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f2_bindings_signals_fdr/f4_cancertype_proper_ctrl'
#     filter_const_file = filter_const_dir + os.sep+ 'expand200k_fdr_LT_0.01_log2FC_LT_1_const_NONoverlapped_{}.bed'.format(cancertype)
#     with open(filter_const_file) as filter_const_inf:
#         filter_const_df = pd.read_csv(filter_const_inf,sep='\t',index_col=3,header=None)
#     return filter_const_df.index





######################################
### CTCF-gene pair correlations
######################################

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
    gene_collections = return_genome_gene_collections()
    dis_df = dis_df[dis_df.columns.intersection(gene_collections)]
    domain_df = domain_df[domain_df.columns.intersection(gene_collections)]
    cor_df = cor_df[cor_df.columns.intersection(gene_collections)]    
    return dis_df,domain_df,cor_df  
#     df = cor_df[domain_df[(dis_df>=a)&(dis_df<=b)]==2]


######################################
### Cancer specific DEG
######################################

def return_cancertype_deg(dataType,fc_thre=0.585,padj_thre=0.01):
    # deg identified by final collection, use for ctcf-gene correlation
    dataTypes = ['AML','BRCA','CRC','LUAD','PRAD','T_ALL','T_ALL_CUTLL1','T_ALL_Jurkat','PD9','PD30']
    
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f6_panCancer_DEG'
    deseq2_file = pardir+os.sep+'/{}/f4_deseq_out_shrink/treated_cancer_vs_ctrl_normal.csv'.format(dataType)
    if dataType=='PD9' or dataType == 'PD30':
        pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f1_TALL_CTCF_binding_vs_GeneExpr/f2_combine_patients/f4_salmon_deseq2_patient/f3_deseq_out/'
        deseq2_file = pardir+os.sep+'treated_{}_vs_ctrl_CD4.csv'.format(dataType)

    with open(deseq2_file) as deseq2_inf:
        deseq2_df = pd.read_csv(deseq2_inf,index_col=0) 
    deseq2_up_genes = deseq2_df[(deseq2_df['log2FoldChange']>fc_thre)&(deseq2_df['padj']<padj_thre)].index.values
    deseq2_down_genes = deseq2_df[(deseq2_df['log2FoldChange']<-1*fc_thre)&(deseq2_df['padj']<padj_thre)].index.values 
    return deseq2_up_genes,deseq2_down_genes




######################################
### specific CTCF intra-domain genes
######################################

# def return_intra_domain_genes(cancertype):
#     # cancertypes = ['T-ALL','Breast_cancer','Colon_cancer','Lung_cancer']
#     intra_domain_gene_list_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/fx_all_feature_analysis_T_ALL/f3_TALL_CTCF_IntraDomain_genes/f5_panCancer_CTCF_intra_domain_genes"
#     gained_gene_file = intra_domain_gene_list_dir+os.sep+'{}_gained_CTCF_intra_domain_genes.csv.txt'.format(cancertype)
#     lost_gene_file = intra_domain_gene_list_dir+os.sep+'{}_lost_CTCF_intra_domain_genes.csv.txt'.format(cancertype)
#     const_gene_file = intra_domain_gene_list_dir+os.sep+'{}_const_CTCF_intra_domain_genes.csv.txt'.format(cancertype)
#     gained_DomainGenes = [i.strip() for i in open(gained_gene_file).readlines()]
#     lost_DomainGenes = [i.strip() for i in open(lost_gene_file).readlines()]
#     const_DomainGenes = [i.strip() for i in open(const_gene_file).readlines()]
#     return gained_DomainGenes,lost_DomainGenes,const_DomainGenes
# 
# 
# def return_notch_occupied_intra_domain_genes(cancertype):
#     # cancertypes = ['T-ALL','Breast_cancer','Colon_cancer','Lung_cancer']
#     intra_domain_gene_list_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/fx_all_feature_analysis_T_ALL/f3_TALL_CTCF_IntraDomain_genes/f6_panCancer_CTCF_notch_occupied_intra_domain_genes"
#     gained_gene_file = intra_domain_gene_list_dir+os.sep+'{}_gained_CTCF_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     lost_gene_file = intra_domain_gene_list_dir+os.sep+'{}_lost_CTCF_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     const_gene_file = intra_domain_gene_list_dir+os.sep+'{}_const_CTCF_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     gained_NotchDomainGenes = [i.strip() for i in open(gained_gene_file).readlines()]
#     lost_NotchDomainGenes = [i.strip() for i in open(lost_gene_file).readlines()]
#     const_NotchDomainGenes = [i.strip() for i in open(const_gene_file).readlines()]
#     return gained_NotchDomainGenes,lost_NotchDomainGenes,const_NotchDomainGenes

# def return_dynamic_notch_occupied_intra_domain_genes(cancertype):
#     # cancertypes = ['T-ALL','Breast_cancer','Colon_cancer','Lung_cancer']
#     intra_domain_gene_list_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/fx_all_feature_analysis_T_ALL/f3_TALL_CTCF_IntraDomain_genes/f6_panCancer_CTCF_notch_occupied_intra_domain_genes"
#     gained_gene_file = intra_domain_gene_list_dir+os.sep+'{}_gained_CTCF_dynamic_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     lost_gene_file = intra_domain_gene_list_dir+os.sep+'{}_lost_CTCF_dynamic_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     const_gene_file = intra_domain_gene_list_dir+os.sep+'{}_const_CTCF_dynamic_notch_occupied_intra_domain_genes.csv.txt'.format(cancertype)
#     gained_NotchDomainGenes = [i.strip() for i in open(gained_gene_file).readlines()]
#     lost_NotchDomainGenes = [i.strip() for i in open(lost_gene_file).readlines()]
#     const_NotchDomainGenes = [i.strip() for i in open(const_gene_file).readlines()]
#     return gained_NotchDomainGenes,lost_NotchDomainGenes,const_NotchDomainGenes

######################################
### GSI/shCTCF/Notch-target genes(DEG)
######################################

def return_logFC_padj(deseq2_file,fcthre,pthre,flag):
    with open(deseq2_file) as inf:
        df = pd.read_csv(inf,index_col=0)
    upgenes = df[(df['log2FoldChange']>fcthre)&(df['pvalue']<pthre)].index.values
    dngenes = df[(df['log2FoldChange']<-1*fcthre)&(df['pvalue']<pthre)].index.values
    print('{}\tfcthre:{}\tpthre:{}\t#upgenes:{}\t#dngenes{}\n'.format(flag,fcthre,pthre,len(upgenes),len(dngenes)))
    return df.index,upgenes,dngenes

# # CUTLL1 shCTCF
def return_CUTLL1_shCTCF_deg(fcthre=0.58,pthre=0.01):
    # return DEG in shCTCF
    #deseq2_CUTLL1_shCTCF = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f3_deseq_out/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
    deseq2_CUTLL1_shCTCF = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f0_processing/salmon_Deseq2/salmon_Deseq2_pca/CUTLL1/f4_deseq_out_shrink/treated_CUTLL1_shCTCF_vs_ctrl_CUTLL1_PIG.csv"
    shCTCF_all,shCTCF_upgenes,shCTCF_dngenes = return_logFC_padj(deseq2_CUTLL1_shCTCF,fcthre,pthre,'shCTCF')
    return shCTCF_all,shCTCF_upgenes,shCTCF_dngenes


# T-ALL cell lines vs. CD4
def return_CUTLL1_vs_CD4_deg(fcthre=0.58,pthre=0.01):
    deseq2_CUTLL1_CD4 = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f2_salmon_deseq2_pca/f3_deseq_out/treated_CUTLL1_vs_ctrl_CD4.csv"
    cutll1_all,cutll1_upgenes,cutll1_dngenes = return_logFC_padj(deseq2_CUTLL1_CD4,fcthre,pthre,'CUTLL1_CD4')
    return cutll1_all,cutll1_upgenes,cutll1_dngenes
# 
def return_Jurkat_vs_CD4_deg(fcthre=0.58,pthre=0.01):
    deseq2_Jurkat_CD4 = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f2_salmon_deseq2_pca/f3_deseq_out/treated_Jurkat_vs_ctrl_CD4.csv"
    jurkat_all,jurkat_upgenes,jurkat_dngenes = return_logFC_padj(deseq2_Jurkat_CD4,fcthre,pthre,'JURKAT_CD4')
    return jurkat_all,jurkat_upgenes,jurkat_dngenes

def return_TALL_vs_CD4_deg(fcthre=0.58,pthre=0.01):
    deseq2_CUTLL1_CD4 = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f2_salmon_deseq2_pca/f3_deseq_out/treated_CUTLL1_vs_ctrl_CD4.csv"
    cutll1_all,cutll1_upgenes,cutll1_dngenes = return_logFC_padj(deseq2_CUTLL1_CD4,fcthre,pthre,'CUTLL1_CD4')
    deseq2_Jurkat_CD4 = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f2_salmon_deseq2_pca/f3_deseq_out/treated_Jurkat_vs_ctrl_CD4.csv"
    jurkat_all,jurkat_upgenes,jurkat_dngenes = return_logFC_padj(deseq2_Jurkat_CD4,fcthre,pthre,'JURKAT_CD4')
    tall_all = np.append(cutll1_all,jurkat_all)
    tall_upgenes = np.append(cutll1_upgenes,jurkat_upgenes)
    tall_dngenes = np.append(cutll1_dngenes,jurkat_dngenes)
    return np.unique(tall_all),np.unique(tall_upgenes),np.unique(tall_dngenes)
 
# NOTCH1 target genes
def return_Notch1_target_genes():
    notch_targets_genelist = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/18_MYC-ChIP_shCTCF-RNA/shCTCF_RNA/f5_shCTCF_DEG_vs_NOTCH1_targets_vs_TALL_DEG/data/notch1_target_dynamic_RP_and_GSIwo_up.txt'
    notch_targets = [i.strip() for i in open(notch_targets_genelist).readlines()]
    return notch_targets


# GSI DEG
def return_CUTLL1_GSI_deg(fcthre=0.58,pthre=0.01,flag='ori'):
    # flag=['ori','trimmed']
    deseq2_CUTLL1_GSI = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_{}/salmon_deseq2_pca_cutll1/f4_deseq_out_shrink/treated_CUTLL1_GSI_vs_ctrl_CUTLL1_DMSO.csv'.format(flag)
    genes_all,genes_upgenes,genes_dngenes = return_logFC_padj(deseq2_CUTLL1_GSI,fcthre,pthre,'CUTLL1_GSI_vs_DMSO')
    return genes_all,genes_upgenes,genes_dngenes
def return_CUTLL1_GSI_wo_deg(fcthre=0.58,pthre=0.01,flag='ori'):
    # flag=['ori','trimmed']
    deseq2_CUTLL1_GSI_wo = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_{}/salmon_deseq2_pca_cutll1/f4_deseq_out_shrink/treated_CUTLL1_GSI_wo16hr_vs_ctrl_CUTLL1_GSI.csv'.format(flag)
    genes_all,genes_upgenes,genes_dngenes = return_logFC_padj(deseq2_CUTLL1_GSI_wo,fcthre,pthre,'CUTLL1_GSI_wo_vs_GSI')
    return genes_all,genes_upgenes,genes_dngenes
# 
# 
# def return_Jurkat_GSI_deg(fcthre=0.58,pthre=0.01,flag='ori'):
#     # flag=['ori','trimmed']
#     deseq2_Jurkat_GSI = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_{}/salmon_deseq2_pca_jurkat/f4_deseq_out_shrink/treated_Jurkat_GSI_vs_ctrl_Jurkat_DMSO.csv'.format(flag)
#     genes_all,genes_upgenes,genes_dngenes = return_logFC_padj(deseq2_Jurkat_GSI,fcthre,pthre,'Jurkat_GSI_vs_DMSO')
#     return genes_all,genes_upgenes,genes_dngenes
# def return_Jurkat_GSI_wo_deg(fcthre=0.58,pthre=0.01,flag='ori'):
#     # flag=['ori','trimmed']
#     deseq2_Jurkat_GSI_wo = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_{}/salmon_deseq2_pca_jurkat/f4_deseq_out_shrink/treated_Jurkat_GSI_wo16hr_vs_ctrl_Jurkat_GSI.csv'.format(flag)
#     genes_all,genes_upgenes,genes_dngenes = return_logFC_padj(deseq2_Jurkat_GSI_wo,fcthre,pthre,'Jurkat_GSI_wo_vs_GSI')
#     return genes_all,genes_upgenes,genes_dngenes




def main():

    pass
    

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
