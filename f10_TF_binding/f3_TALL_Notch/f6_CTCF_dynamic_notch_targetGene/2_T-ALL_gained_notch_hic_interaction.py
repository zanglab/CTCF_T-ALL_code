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



def process_df_get_hic_interaction(df,interaction_file,outfile_name):

    #print(df,df.columns);exit()
    resolution = 5000
    with open(interaction_file) as interaction_inf:
        interaction_df = pd.read_csv(interaction_inf,sep='\t',index_col=0)
    print(outfile_name)
    outf = open(outfile_name,'w')
    outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('id','chr','ctcf_pos','notch_pos','ctcf_notch_interaction','intra_domain_ctcf_all_interactions'))
    for id_index in df.index:
        chr = df.loc[id_index,'chr']
        domain_left = df.loc[id_index,'domain_100k_1M_left']
        domain_right = df.loc[id_index,'domain_100k_1M_right']
        pos1 = df.loc[id_index,'mid_position']
        overlap_info = df.loc[id_index,'if_intra_domain_dynamic_notch']
        #if overlap_info==1 and id_index in interaction_df.index:
        if overlap_info==1:
            # get all 2m interaction score of this binding site
            interactions_2m_list = interaction_df.loc[id_index]
            # adjust position by resolution
            adj_pos1 = pos1//resolution*resolution
            adj_domain_left = domain_left//resolution*resolution
            adj_domain_right = domain_right//resolution*resolution
            # here is the pos in interactions_2m_list
            domain_left_list_pos = adj_domain_left - adj_pos1
            domain_right_list_pos = adj_domain_right - adj_pos1
            a = [str(int(i)) for i in np.arange(domain_left_list_pos,domain_right_list_pos+resolution,resolution)]
            bg_scores = interactions_2m_list.loc[a].values
            bg_scores = ','.join([str(i) for i in bg_scores])
            #print('id,ctcf_pos:\t',id_index,pos1,adj_pos1)
            #print('domain-left:\t',domain_left,adj_domain_left)
            #print('domain-right:\t',domain_right,adj_domain_right)
            #print(bg_scores)
            
            notch_poses = [i for i in map(int,df.loc[id_index,'dynamic_notch_center'].split(','))]
            compr_score=0
            notch_pos = 0
            for pos2 in notch_poses:
                # adj and list pos for notch binding
                adj_pos2 = pos2//resolution*resolution
                compr_list_pos = int(adj_pos2 - adj_pos1)#;print(pos2,adj_pos2,compr_list_pos)
                tmp_score = interactions_2m_list.loc[str(compr_list_pos)];print(id_index,pos2,tmp_score)
                outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(id_index,chr,pos1,pos2,tmp_score,bg_scores))
#                 if tmp_score>=compr_score:
#                     compr_score = tmp_score
#                     notch_pos = pos2  
#             #print(compr_score,notch_pos)   
#             outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(id_index,chr,pos1,notch_pos,compr_score,bg_scores))
    outf.close()



def main():

    outdir = 'f2_CTCF_Notch1_hic_interaction'
    os.makedirs(outdir,exist_ok=True) 
    a6010_raw_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f2_union_binding_view2M_bothside_interaction/A6010_raw_res5000_union.csv'
    cutll1_raw_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f2_union_binding_view2M_bothside_interaction//Jurkat_raw_res5000_union.csv'
    cutll1_raw_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f4_HiC_TADs/f1_union_binding_processed/f2_union_binding_view2M_bothside_interaction//Cutll1_raw_res5000_union.csv'

    union_feature_df = CTCF_TALL_modules_new.return_cancer_specific_combined_features('T-ALL')
    gained_df,lost_df = CTCF_TALL_modules_new.return_cancer_specific_binding('T-ALL')

    #
    basename = 'T-ALL_gained'
    df_tmp = union_feature_df.loc[gained_df.index]
    process_df_get_hic_interaction(df_tmp,a6010_raw_file,outdir+os.sep+'{}_ctcf_notch_interaction_vs_bg_in_{}.csv'.format(basename,'a6010_raw'))
    process_df_get_hic_interaction(df_tmp,cutll1_raw_file,outdir+os.sep+'{}_ctcf_notch_interaction_vs_bg_in_{}.csv'.format(basename,'cutll1_raw'))

    #
    basename = 'T-ALL_lost'
    df_tmp = union_feature_df.loc[lost_df.index]
    process_df_get_hic_interaction(df_tmp,a6010_raw_file,outdir+os.sep+'{}_ctcf_notch_interaction_vs_bg_in_{}.csv'.format(basename,'a6010_raw'))
    process_df_get_hic_interaction(df_tmp,cutll1_raw_file,outdir+os.sep+'{}_ctcf_notch_interaction_vs_bg_in_{}.csv'.format(basename,'cutll1_raw'))
    



 
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

