# <2019-01-26>
'''
@ authors: Zhenjia Wang, Chongzhi Zang

For each region in bed.file, keep all UDHS regions overlapped with this region
For each TF, keep the binding info on all overlapped UDHS regions for each region
Then for each TF, compare the overlap enrichment .bed files compare to all UDHS regions
'''


import sys,argparse
import os,glob,re
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
import CTCF_TALL_modules
import scipy
from scipy import stats
import overlapped_udhs_id
import stats_ctrl_UDHS_hg38


def return_udhs_overlapped_ids(df_dic,motif_file):
    # read df from each TF-UDHS file
    tfname = os.path.basename(motif_file).split('_')[0] ;print(tfname)
    tf = open(motif_file, 'rb') 
    lines = tf.raw.read()
    tf_dic = {}
    for key in df_dic.keys():
        match = [lines[2*position-2]-ord('0') for position in df_dic[key]]
        tf_dic[key]=sum(match)
    tf.close()
    #print(tfname,tf_dic);exit()
    return tfname,tf_dic
    


def bed_file_udhs_overlapped_id_group_df(udhs_overlapped):
    df = pd.DataFrame()
    df_dic = {}
    # for each region in the bed file, record overlapped UDHS. One region may overlap with multiple UDHS regions
    for chr in udhs_overlapped.keys():
        for overlapped_region in udhs_overlapped[chr]:
            start,end = overlapped_region[0],overlapped_region[1]
            overlapped_ids = overlapped_region[2].split('\t')[-1]
            ## keep the dic info to get the 0/1 overlap info from each TF file
            overlapped_ids_list = [int(i) for i in overlapped_region[2].split('\t')[-1].split(',')]
            #print(overlapped_ids)
            index_id='{}_{}_{}'.format(chr,start,end)
            df.loc[index_id,'overlapped_udhs'] = overlapped_ids
            df_dic[index_id] = overlapped_ids_list
    #print(df_dic);exit()
    df = df.fillna(0)
    return df,df_dic



def main(treat,control,outdir):

    #outdir = 'f1_motif_overlap_csv'
    os.makedirs(outdir,exist_ok=True) 
    udhs_file = '/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
    udhs_motif_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2017/tf_merged/hg38'
    #udhs_motif_dir="/nv/vol190/zanglab/zw5j/data/ChIPseq_overlap_2016_new/TF_DATA_455"
    udhs_motif_files = glob.glob(udhs_motif_dir+'/*.txt')
            
    udhs_overlapped,udhs_nonoverlapped = overlapped_udhs_id.return_udhs_overlapped_ids(treat,udhs_file,'hg38')
    out_df,df_dic = bed_file_udhs_overlapped_id_group_df(udhs_overlapped)
    #print(udhs_overlapped,'\n\n',udhs_nonoverlapped);exit()
    basename = os.path.basename(treat).split('.bed')[0]#;print(basename)
    for motif_file in sorted(udhs_motif_files):
        tfname,tf_dic = return_udhs_overlapped_ids(df_dic,motif_file)
        out_df = pd.concat([out_df,pd.Series(tf_dic).to_frame(tfname)],axis=1)        
    out_df.to_csv(outdir+os.sep+'{}_udhs_chip_peak_overlapped.csv'.format(basename))
    stats_ctrl_UDHS_hg38.stats_test(out_df,basename,outdir)

    if control:
        udhs_overlapped,udhs_nonoverlapped = overlapped_udhs_id.return_udhs_overlapped_ids(control,udhs_file,'hg38')
        control_out_df,control_df_dic = bed_file_udhs_overlapped_id_group_df(udhs_overlapped)
        control_basename = os.path.basename(control).split('.bed')[0]#;print(basename)
        for motif_file in sorted(udhs_motif_files):
            tfname,tf_dic = return_udhs_overlapped_ids(control_df_dic,motif_file)
            control_out_df = pd.concat([control_out_df,pd.Series(tf_dic).to_frame(tfname)],axis=1)        
        control_out_df.to_csv(outdir+os.sep+'{}_udhs_chip_peak_overlapped.csv'.format(control_basename))
        
        stats_ctrl_UDHS_hg38.stats_test_with_control(out_df,control_out_df,basename,outdir)




 
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--treat', action = 'store', type = str,dest = 'treat', help = 'input bed file of treat', metavar = '<file>')
    parser.add_argument('-c','--control', action = 'store', type = str,dest = 'control', help = 'input bed file of control', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.treat,args.control,args.outdir)

