import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
# from GenomeData import *
from operator import itemgetter
#def expand_region(summitlist):
import CTCF_TALL_modules_new 
# import association_with_regions
from collections import Counter


def main():               
    
    outdir='f2_tmp'
    os.makedirs(outdir,exist_ok=True)
    
#     union_df = CTCF_TALL_modules_new.return_union_binding_df()
    motif_df = CTCF_TALL_modules_new.return_union_binding_with_motif()
    occupancy_df = CTCF_TALL_modules_new.return_occupancy_df()
#     gene_annotation_file = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f5_add_gene_annotation/union_binding_with_gene_annotation.csv'
#     gene_annotation_df = pd.read_csv(gene_annotation_file,index_col=3)
    
    # get the combined df and mid positions
    combined_df = pd.concat([motif_df,occupancy_df['occupancy_score']],axis=1)
    combined_df.insert(3,'id',combined_df.index)
#     combined_df['annotation'] = gene_annotation_df['annotation']
    
    motif_mid = 0.5*(combined_df['motif_start'].replace("N",0).astype(int)+combined_df['motif_end'].replace("N",0).astype(int))
    binding_mid = 0.5*(combined_df['start']+combined_df['end'])
    combined_df['mid_position'] = pd.concat([binding_mid[motif_mid==0],motif_mid[motif_mid!=0]]).astype(int)

    # save the combined data
    combined_df.to_csv('{}/union_binding_tmp.csv'.format(outdir),index=False)
    combined_df.iloc[:,:4].to_csv('{}/union_binding_tmp.bed'.format(outdir),index=False,sep='\t',header=None)
    
#     combined_df_motif = combined_df[combined_df['motif_strand']!='N']
#     combined_df_motif.to_csv('{}/union_binding_with_motif.csv'.format(outdir),index=False)
#     combined_df_motif.iloc[:,:4].to_csv('{}/union_binding_with_motif.bed'.format(outdir),index=False,sep='\t',header=None)
#     
#     # high quality bindings (occupancy score ≥3)
#     occupancy3_filtered = combined_df[combined_df['occupancy_score']>=3]
#     occupancy3_filtered.to_csv('{}/union_binding_occupancy_score_GT3.csv'.format(outdir),index=False)
#     occupancy3_filtered.iloc[:,:4].to_csv('{}/union_binding_occupancy_score_GT3.bed'.format(outdir),index=False,sep='\t',header=None)
#     
#     occupancy3_filtered_motif = occupancy3_filtered[occupancy3_filtered['motif_strand']!='N']
#     occupancy3_filtered_motif.to_csv('{}/union_binding_occupancy_score_GT3_with_motif.csv'.format(outdir),index=False)
#     occupancy3_filtered_motif.iloc[:,:4].to_csv('{}/union_binding_occupancy_score_GT3_with_motif.bed'.format(outdir),index=False,sep='\t',header=None)
#     
#     # constitutive bindings
#     constitutive_thre = 0.90
#     constitutive_df = combined_df[combined_df['occupancy_score']>=771*constitutive_thre]
#     constitutive_df.to_csv('{}/union_binding_constitutive_{}.csv'.format(outdir,constitutive_thre),index=False)
#     constitutive_df.iloc[:,:4].to_csv('{}/union_binding_constitutive_{}.bed'.format(outdir,constitutive_thre),index=False,sep='\t',header=None)
#     
#     constitutive_thre = 0.85
#     constitutive_df = combined_df[combined_df['occupancy_score']>=771*constitutive_thre]
#     constitutive_df.to_csv('{}/union_binding_constitutive_{}.csv'.format(outdir,constitutive_thre),index=False)
#     constitutive_df.iloc[:,:4].to_csv('{}/union_binding_constitutive_{}.bed'.format(outdir,constitutive_thre),index=False,sep='\t',header=None)
#     
#     constitutive_thre = 0.80
#     constitutive_df = combined_df[combined_df['occupancy_score']>=771*constitutive_thre]
#     constitutive_df.to_csv('{}/union_binding_constitutive_{}.csv'.format(outdir,constitutive_thre),index=False)
#     constitutive_df.iloc[:,:4].to_csv('{}/union_binding_constitutive_{}.bed'.format(outdir,constitutive_thre),index=False,sep='\t',header=None)
#     
#     constitutive_thre = 0.75
#     constitutive_df = combined_df[combined_df['occupancy_score']>=771*constitutive_thre]
#     constitutive_df.to_csv('{}/union_binding_constitutive_{}.csv'.format(outdir,constitutive_thre),index=False)
#     constitutive_df.iloc[:,:4].to_csv('{}/union_binding_constitutive_{}.bed'.format(outdir,constitutive_thre),index=False,sep='\t',header=None)
#     
    
    
    
        
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()