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
import association_with_genes
import association_with_regions
import re,bisect
import CTCF_TALL_modules_new
import scipy
from scipy import stats


cancer_type_hic_filename = {'T-ALL':['Jurkat','A6010'],\
                                'T-ALL-P1':['PD9','A6010'],\
                                'T-ALL-P2':['PD31','A6010'],\
                                'CRC':['HCT116','trans_colon1'],\
                                'PRAD':['LNCaP','PrEC'],\
                                'PRAD_TissueAdded':['LNCaP','PrEC']}
    

def return_hic_file(cancertype,normalization,resolution):
    union_hic_file_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f2_union_binding_view2M_bothside_interaction'
    cancer_union_hic_file = union_hic_file_dir+os.sep+'{}_{}_res{}_union.csv'.format(cancer_type_hic_filename[cancertype][0],normalization,resolution)
    normal_union_hic_file = union_hic_file_dir+os.sep+'{}_{}_res{}_union.csv'.format(cancer_type_hic_filename[cancertype][1],normalization,resolution)
    return cancer_union_hic_file,normal_union_hic_file
    

def process_df_get_hic_interaction(df,interaction_file,outfile_name):
    # == get df using CTCF id
    with open(interaction_file) as interaction_inf:
        interaction_df = pd.read_csv(interaction_inf,sep='\t',index_col=0)
    df = interaction_df.loc[df.index]
    df.to_csv(outfile_name)


    
def return_cancertype_abstraction(cancertype,normalization,resolution):
    outdir = 'f1_gained_lost_binding_interactions'
    os.makedirs(outdir,exist_ok=True) 
    cancer_hic_file,normal_hic_file = return_hic_file(cancertype,normalization,resolution)
    print(os.path.basename(normal_hic_file))
        
    if cancertype in ['T-ALL-P1','T-ALL-P2']:
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')
    else:
        gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding(cancertype)
    #
    basename = '{}_gained'.format(cancertype)
    process_df_get_hic_interaction(gained_df,cancer_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][0],normalization,resolution))
    process_df_get_hic_interaction(gained_df,normal_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][1],normalization,resolution))

    #
    basename = '{}_lost'.format(cancertype)
    process_df_get_hic_interaction(lost_df,cancer_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][0],normalization,resolution))
    process_df_get_hic_interaction(lost_df,normal_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][1],normalization,resolution))
    
    #
    constitutive_df = CTCF_TALL_modules_new.return_constitutive_df()#;print(constitutive_df);exit()
    basename = '{}_ctrl'.format(cancertype)
    process_df_get_hic_interaction(constitutive_df,cancer_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][0],normalization,resolution))
    process_df_get_hic_interaction(constitutive_df,normal_hic_file,outdir+os.sep+'{}_ctcf_interactions_in_{}_{}_res{}.csv'.format(basename,cancer_type_hic_filename[cancertype][1],normalization,resolution))
    


def main():

    cancertypes=['T-ALL','T-ALL-P1','T-ALL-P2','CRC','PRAD','PRAD_TissueAdded']
    cancertypes=['T-ALL-P1','T-ALL-P2','CRC','PRAD','PRAD_TissueAdded']
    for cancertype in cancertypes:
        for normalization in  ['raw']:
            for resolution in [5000]:
                return_cancertype_abstraction(cancertype,normalization,resolution)




 
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

