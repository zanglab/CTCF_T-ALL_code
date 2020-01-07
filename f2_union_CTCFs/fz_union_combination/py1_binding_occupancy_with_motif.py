#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 22:00:55 2017
@author: zw5j
"""

import os,sys,re,argparse,bisect
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.6)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
import CTCF_TALL_modules_new


def return_motif_only_binding_df():
    motif_bed = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f4_motif_on_union_CTCFs/f1_union_binding_motif/union_CTCF_with_motif_with_motif_only.bed'  
    with open(motif_bed) as motif_inf:
        motif_df = pd.read_csv(motif_inf,sep='\t',index_col=3)
    return motif_df

def return_occurrence_binding_df():
    input_bed = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/f2_occupancy_on_union_CTCFs/f2_occupancy_score/union_CTCF_occupancy_score.csv'  
    with open(input_bed) as input_inf:
        input_df = pd.read_csv(input_inf,index_col=3,sep='\t')
    return input_df
    

def main():

    collection_df = CTCF_TALL_modules_new.return_collection_df()
    occurrence_df = return_occurrence_binding_df()
    motif_df = return_motif_only_binding_df()
    
    out_df = pd.DataFrame()
    for i in np.arange(1,collection_df.shape[0]+1):
        occurred_ids = occurrence_df[occurrence_df['occupancy_score']==i].index
        occurred_ids_with_motif = occurred_ids.intersection(motif_df.index)
        out_df.loc[i,'#bindings']=len(occurred_ids)
        out_df.loc[i,'#bindings_with_motif']=len(occurred_ids_with_motif)
        print(i)
    out_df['%motif']= np.around(out_df['#bindings_with_motif']/out_df['#bindings'],2) 
    out_df.to_csv('binding_occupancy_with_motif.csv')
    




if __name__ == '__main__':
  
    main()

