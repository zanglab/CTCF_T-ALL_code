'''
select within domain interactions between jurkat and cd4 then cal logFC
'''

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
# import association_with_genes
# import association_with_regions
import re,bisect
# import CTCF_TALL_modules
import scipy
from scipy import stats
import MA_regression


cancer_type_hic_filename = {'T-ALL':['Jurkat','A6010'],\
                                'T-ALL-P1':['PD9','A6010'],\
                                'T-ALL-P2':['PD31','A6010'],\
                                'CRC':['HCT116','trans_colon1'],\
                                'PRAD':['LNCaP','PrEC'],\
                                'PRAD_TissueAdded':['LNCaP','PrEC']}

## ====


def hic_cor_high_logFC_into_bed(df,cancer_interaction_file,normal_interaction_file,outdir,cancertype,data_type):
    #print(os.path.isfile(interaction_file))#;exit()
    resolution = 5000
    with open(cancer_interaction_file) as cancer_interaction_inf:
        cancer_df = pd.read_csv(cancer_interaction_inf,index_col=0)
    with open(normal_interaction_file) as normal_interaction_inf:
        normal_df = pd.read_csv(normal_interaction_inf,index_col=0)
    assert all(cancer_df.index == normal_df.index)
    df = df.loc[cancer_df.index]
    df = df.dropna()
    
    combined_df = pd.DataFrame()
    for id_index in df.index:
        chr = df.loc[id_index,'chr']
        domain_left = df.loc[id_index,'domain_100k_1M_left']
        domain_right = df.loc[id_index,'domain_100k_1M_right']
        pos1 = df.loc[id_index,'middle']
        adj_pos1 = pos1//resolution*resolution
        adj_domain_left = domain_left//resolution*resolution
        adj_domain_right = domain_right//resolution*resolution
        # here is the pos in interactions_2m_list
        domain_left_list_pos = adj_domain_left - adj_pos1
        domain_right_list_pos = adj_domain_right - adj_pos1
        a = [str(int(i)) for i in np.arange(domain_left_list_pos,domain_right_list_pos+resolution,resolution)]
#         print(cancertype,len(a),type(domain_left_list_pos),domain_right_list_pos)
        #a = [str(int(i)) for i in np.arange(-500000,500000+resolution,resolution)]
        cancer_domain_list = cancer_df.loc[id_index].loc[a]#.value
#         print(cancer_domain_list);exit()
        normal_domain_list = normal_df.loc[id_index].loc[a]#.values
        compr_df = pd.concat([cancer_domain_list,normal_domain_list],axis=1)
        compr_columns = cancer_type_hic_filename[cancertype]
        compr_df.columns = compr_columns
        compr_df.index = '{}_{}_'.format(chr,adj_pos1)+compr_df.index
        # combine all info from all bindings
        combined_df = pd.concat([compr_df,combined_df])
        
    combined_df = MA_regression.ma_fit(combined_df,compr_columns)
    #print(combined_df);exit()
    outfile_name = outdir+os.sep+'{}_{}_logFC.bed'.format(cancertype,data_type)
    combined_df.to_csv(outfile_name)



def return_domain_df():
    domain_file = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo.csv"
    domain_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f6_gene_expr/f4_gene_CTCF_matrices/f2_gene_CTCF_domain/all_CTCF_domainInfo_GT100K_LT1M_EachSide.csv'
    with open(domain_file) as domain_inf:
        domain_df = pd.read_csv(domain_inf,sep=',',index_col=2)
    return domain_df

def return_cancer_specific_binding(cancertype):
    # == cancer specific gained/lost bindings ====
    pardir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f3_cancer_specific_CTCFs/f3_cancer_specific_gained_lost/f1_cancer_specific_binding'
    cancertypes=['T-ALL','BRCA','CRC','LUAD','AML','PRAD','PRAD_TissueAdded']
    gained_file = pardir+os.sep+'{}_gained.csv'.format(cancertype)
    lost_file = pardir+os.sep+'{}_lost.csv'.format(cancertype)
    with open(gained_file) as gained_inf, open(lost_file) as lost_inf:
        gained_df = pd.read_csv(gained_inf,index_col=0)
        lost_df = pd.read_csv(lost_inf,index_col=0)
    return gained_df,lost_df


def main():

    outdir = 'f1_high_cor_high_logFC_regions'
    os.makedirs(outdir,exist_ok=True) 

    domain_df = return_domain_df()
    hic_interaction_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f2_gained_lost_hic_interaction/f1_gained_lost_binding_interactions"
        
    for cancertype in ['T-ALL','CRC']:
        for data_type in ['gained','lost']:
            cancer_interaction_file = hic_interaction_dir+os.sep+'{}_{}_ctcf_interactions_in_{}_raw_res5000.csv'.format(cancertype,data_type,cancer_type_hic_filename[cancertype][0])
            normal_interaction_file = hic_interaction_dir+os.sep+'{}_{}_ctcf_interactions_in_{}_raw_res5000.csv'.format(cancertype,data_type,cancer_type_hic_filename[cancertype][1])
#             print(os.path.basename(cancer_interaction_file),os.path.basename(normal_interaction_file))
            hic_cor_high_logFC_into_bed(domain_df,cancer_interaction_file,normal_interaction_file,outdir,cancertype,data_type)



 
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

